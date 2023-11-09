=head

BEFORE USAGE, RUN:

module add python-cbrg/202111
module add R-base/4.1.0
module add R-cbrg/202111 # R-base/4.1.0 gsl/2.6 hdf5/1.10.7
conda activate kbpython

=cut

use strict;
use warnings;

# REQUIREMENTS
my $root       = '/project/GorielyLab2021/sbush/ssc_atlas';
my $species    = 'Homo_sapiens'; # 'Bos_grunniens'; # 'Macaca_mulatta'; # 'Papio_anubis'; # 'Mus_musculus'; # 'Gallus_gallus'; # 'Ailuropoda_melanoleuca'; # 'Bubalus_bubalis'; # 'Macaca_fascicularis'; # 'Ovis_aries'; # 'Sus_scrofa'; # 'Rattus_norvegicus';
my $fq_dir     = "/project/GorielyLab2021/sbush/ssc_atlas/fq/$species"; # from 1.download_fqs.pl
my $t2g        = "$root/indexes/$species/t2g2.txt"; # from 3.update_t2g_to_replace_null_names.pl
my $index      = "$root/indexes/$species/index.idx"; # from 2.make_kb_index.sh

# input files for human RNA-seq velocity analysis
#my $t2g           = "$root/indexes/velocity/t2g2.txt"; # from 3.update_t2g_to_replace_null_names.pl
#my $index         = "$root/indexes/velocity/index.idx"; # from 2.make_kb_index.sh
#my $intron_fa     = "$root/indexes/velocity/intron.fa"; # from 2.make_kb_index.sh
#my $spliced_t2c	  = "$root/indexes/velocity/cdna_t2c.txt"; # from 2.make_kb_index.sh
#my $unspliced_t2c = "$root/indexes/velocity/intron_t2c.txt"; # from 2.make_kb_index.sh

my $metadata   = "$root/prerequisites/metadata.$species.tsv"; # manually created
my $make_plots = "$root/make_plots_SCTransform.$species.R"; # manually created; see https://www.kallistobus.tools/tutorials/kb_building_atlas/r/kb_analysis_0_r. NOTE THAT THIS R SCRIPT CONTAINS A HARD-CODED PATH TO $t2g
my $fatal      = 0;
if (!(-e($make_plots))) { $fatal++; print "ERROR: cannot find $make_plots\n"; }
if (!(-e($metadata)))   { $fatal++; print "ERROR: cannot find $metadata\n";   }
if (!(-d($fq_dir)))     { $fatal++; print "ERROR: cannot find $fq_dir\n";     }
if (!(-e($index))) 	    { $fatal++; print "ERROR: cannot find $index\n"; 	  }
if (!(-e($t2g))) 	    { $fatal++; print "ERROR: cannot find $t2g\n"; 	      }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 15;

# OUTPUT
my $out_dir = "$root/kb";
if (!(-d($out_dir))) 			{ mkdir $out_dir 			or die $!; }
if (!(-d("$out_dir/$species"))) { mkdir "$out_dir/$species" or die $!; }
my $out_sh = "$root/run_kb.$species.sh";
open(SH,'>',$out_sh) or die $!;
print SH "#!/bin/bash\n";
print SH "mkdir $out_dir/$species\n" unless (-d("$out_dir/$species"));
print SH "cd $out_dir/$species\n";

# WHAT SAMPLE IDs ARE WE GOING TO RUN?
my %samples = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $study = $line[0]; my $sample_id = $line[1]; my $layout = $line[3]; my $age_category = $line[4]; my $techno = $line[6]; my $phenotype = $line[10];
	  $samples{$sample_id}{layout} = $layout;
	  $samples{$sample_id}{techno} = $techno;
	}
close(IN) or die $!;

# FOR EACH SAMPLE...
my @samples = ();
while((my $sample_id,my $irrel)=each(%samples))
	{ push(@samples,$sample_id); }
my @sorted_samples = sort {$a cmp $b} @samples;
my $num = 0; my $total = @sorted_samples;
for(my $x=0;$x<@sorted_samples;$x++)
	{ $num++;
	  print "$num of $total\n";
	  my $sample_id = $sorted_samples[$x];
	  my $techno    = $samples{$sample_id}{techno};
	  my $layout    = $samples{$sample_id}{layout};
	  if ($techno eq 'GEXSCOPE') { $techno = 'DROPSEQ'; } # i.e. the Chen 2020 samples, described in this paper: https://academic.oup.com/nar/article/47/3/e16/5193346. This states that "transcriptome alignment including barcode/UMI identification and collapsing were performed as described in Dropseq method" - so for the purpose of quantification with KB, GEXSCOPE = DROPSEQ

	  next if ($layout ne 'PAIRED'); # CHECKPOINT: restrict analysis to paired-end samples only
	  if (($techno !~ /^10X/) and ($techno !~ /DROPSEQ/) and ($techno !~ /INDROPS/))
		{ print "WARNING: skipping $sample_id as we do not support the $techno platform\n"; }
	  next unless (($techno =~ /^10X/) or ($techno =~ /DROPSEQ/) or ($techno =~ /INDROPS/)); # CHECKPOINT: restrict analysis to Chromium, DropSeq or InDrops data only

	  # CHECK WHETHER WE HAVE THE RAW DATA
	  my $fq1 = ''; my $fq2 = ''; my $fq = ''; my $files_available = 0;
	  if ($layout eq 'PAIRED')
		{ $fq1 = "$fq_dir/$sample_id/$sample_id.1.fq.gz";
		  $fq2 = "$fq_dir/$sample_id/$sample_id.2.fq.gz";
		  if ( (-e($fq1)) and (-e($fq2)) ) { $files_available++; }
		}
	  elsif ($layout eq 'SINGLE')
		{ $fq = "$fq_dir/$sample_id/$sample_id.fq.gz";
		  if (-e($fq)) { $files_available++; }
		}
	  if ($files_available == 0) { print "WARNING: skipping $sample_id as we do not have the fqs\n"; }
	  next if ($files_available == 0);
	  
	  # SKIP IF WE'VE *COMPLETELY* PROCESSED THIS SAMPLE ALREADY. IF WE HAVE, WE'LL HAVE PRODUCED AN OUTPUT DIRECTORY CONTAINING MULTIPLE FILES BUT TWO KEY ONES IN PARTICULAR: (1) A COUNT MATRIX: counts_unfiltered/cells_x_genes.mtx, AND (2) THE SAME MATRIX REPRESENTED AS AN ANNDATA OBJECT (FOR LATER INTEGRATION): counts_unfiltered/adata.h5ad
	  if ( (-e("$out_dir/$species/$sample_id/counts_filtered/cells_x_genes.mtx")) and (-e("$out_dir/$species/$sample_id/counts_filtered/adata.h5ad")) and (-e("$out_dir/$species/$sample_id/R_SCT/$sample_id.rds")) )
		{ print "skipping $sample_id as we've completely processed this sample already\n"; }
	  next if ( (-e("$out_dir/$species/$sample_id/counts_filtered/cells_x_genes.mtx")) and (-e("$out_dir/$species/$sample_id/counts_filtered/adata.h5ad")) and (-e("$out_dir/$species/$sample_id/R_SCT/$sample_id.rds")) );
	  
	  # HAVE WE *INCOMPLETELY* PROCESS THIS SAMPLE?
	  # if we've already produced an output directory and yet haven't produced the cell x gene matrix and h5ad (see above checkpoint), then something has gone wrong - kb has not run to completion. In that case, we shall delete this directory and start again.
	  if ( (-d("$out_dir/$species/$sample_id")) and ( (!(-e("$out_dir/$species/$sample_id/counts_filtered/cells_x_genes.mtx"))) or (!(-e("$out_dir/$species/$sample_id/counts_filtered/adata.h5ad"))) ) )
		{ print SH "rm -r $out_dir/$species/$sample_id\n"; }
	  
	  # RUN KB
	  if ( (!(-e("$out_dir/$species/$sample_id/counts_filtered/cells_x_genes.mtx"))) or (!(-e("$out_dir/$species/$sample_id/counts_filtered/adata.h5ad"))) )
		{ if ($layout eq 'PAIRED')
			{ print SH "kb count -i $index -g $t2g -x $techno -o $sample_id --filter bustools --h5ad -t $num_procs $fq1 $fq2\n"; }
		  elsif ($layout eq 'SINGLE')
			{ print SH "kb count -i $index -g $t2g -x $techno -o $sample_id --filter bustools --h5ad -t $num_procs $fq\n"; }
		  print SH "rm $out_dir/$species/$sample_id/output.bus $out_dir/$species/$sample_id/output.filtered.bus $out_dir/$species/$sample_id/output.unfiltered.bus\n";
		}
	  
	  # RUN KB (FOR VELOCITY ANALYSIS)
#	  if ( (!(-e("$out_dir/$species/$sample_id/counts_filtered/adata.loom"))) or (!(-e("$out_dir/$species/$sample_id/counts_unfiltered/adata.loom"))) )
#		{ if ($layout eq 'PAIRED')
#			{ print SH "kb count -i $index -g $t2g -x $techno -o $sample_id -c1 $spliced_t2c -c2 $unspliced_t2c --workflow lamanno --filter bustools --loom -t $num_procs $fq1 $fq2\n"; }
#		  elsif ($layout eq 'SINGLE')
#			{ print SH "kb count -i $index -g $t2g -x $techno -o $sample_id -c1 $spliced_t2c -c2 $unspliced_t2c --workflow lamanno --filter bustools --loom -t $num_procs $fq\n"; }
#		  print SH "rm $out_dir/$species/$sample_id/matrix.ec\n";
#		  print SH "rm $out_dir/$species/$sample_id/output.bus $out_dir/$species/$sample_id/output.filtered.bus $out_dir/$species/$sample_id/output.unfiltered.bus\n";
#		  print SH "rm $out_dir/$species/$sample_id/spliced.filtered.bus $out_dir/$species/$sample_id/spliced.unfiltered.bus\n";
#		  print SH "rm $out_dir/$species/$sample_id/unspliced.filtered.bus $out_dir/$species/$sample_id/unspliced.unfiltered.bus\n";
#		}

	  # CREATE TABLES AND FIGURE USING R
	  # NOTE THAT THE R SCRIPT CONTAINS A HARD-CODED PATH TO $t2g
	  if (!(-d("$out_dir/$species/$sample_id/R_SCT"))) { print SH "mkdir $out_dir/$species/$sample_id/R_SCT\n"; }
	  if (!(-e("$out_dir/$species/$sample_id/R_SCT/$sample_id.rds")))
		{ print SH "cd $out_dir/$species/$sample_id\n";
		  print SH "Rscript $make_plots\n";
		  print SH "mv $out_dir/$species/$sample_id/R_SCT/testis.rds $out_dir/$species/$sample_id/R_SCT/$sample_id.rds\n";
		  print SH "cd $out_dir/$species\n";
		}
	}
close(SH) or die $!;
exit 1;