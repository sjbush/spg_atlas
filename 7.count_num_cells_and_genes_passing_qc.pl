# PURPOSE: count the number of cells and genes per sample that pass QC and so are eligible for inclusion in the whole-testes atlas

use strict;
use warnings;
use Acme::Tools qw(avg median);
use POSIX qw(ceil);

# REQUIREMENTS
my $root     		= '/project/GorielyLab2021/sbush/ssc_atlas';
my $species  		= 'Homo_sapiens';
my $in_dir   		= "$root/kb/$species"; # from 4.run_kb.pl
my $atlas_metadata  = "$root/human_adult.metadata.txt"; # from 6.refine_clustering_of_whole_testes_atlas.R
my $sample_metadata = "$root/prequisites/metadata.$species.tsv"; # manually created
my $fatal    		= 0;
if (!(-d($in_dir)))   		 { $fatal++; print "ERROR: cannot find $in_dir\n";   		}
if (!(-e($atlas_metadata)))  { $fatal++; print "ERROR: cannot find $atlas_metadata\n";  }
if (!(-e($sample_metadata))) { $fatal++; print "ERROR: cannot find $sample_metadata\n"; }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file = "$root/human_adult.no_of_cells_and_genes_per_sample_passing_qc.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "Sample ID\tPublication\tOriginal total no. of cells (unfiltered by KB)\tNo. of cells passing KB internal filters (this number used as basis of subsequent percentages)\t";
print OUT "No. of cells with <= 1000 genes\t% of cells with <= 1000 genes\tNo. of cells with >= 10,000 genes\t% of cells with >= 10,000 genes\t";
print OUT "No. of cells with <= 2000 UMIs\t% of cells with <= 2000 UMIs\tNo. of cells with >= 50,000 UMIs\t% of cells with >= 50,000 UMIs\t";
print OUT "No. of cells with >= 5% MT UMIs\t% of cells with >= 5% MT UMIs\t";
print OUT "No. of cells passing all QC filters (10,000 > genes > 1000 and 50,000 > UMIs > 2000 and < 5% MT)\t% of cells passing all QC filters (10,000 > genes > 1000 and 50,000 > UMIs > 2000 and < 5% MT)\t";
print OUT "No. of genes in >= 3 cells passing all QC filters\t";
print OUT "Minimum genes/cell in integrated dataset\tMean genes/cell in integrated dataset\tMedian genes/cell in integrated dataset\tMaximum genes/cell in integrated dataset\n";

# WHAT SAMPLE IDs ARE WE GOING TO QUANTIFY?
my %publications = ();
open(IN,$sample_metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $study = $line[0]; my $sample_id = $line[1]; my $layout = $line[3]; my $age_category = $line[4]; my $techno = $line[6]; my $phenotype = $line[10];
	  $publications{$sample_id} = $study;
	}
close(IN) or die $!;

# *IN THE INTEGRATED DATASET*, HOW MANY CELLS (i.e. BARCODES) ARE IN EACH SAMPLE, AND HOW MANY GENES ARE REPRESENTED?
my %genes_per_sample = ();
open(IN,$atlas_metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $num_genes = $line[1]; my $barcode = $line[8]; my $sample_id = $line[9];
	  $genes_per_sample{$sample_id}{$barcode} = $num_genes;
	}
close(IN) or die $!;

opendir(DIR,$in_dir) or die $!;
my @sample_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
foreach my $sample_id (@sorted_sample_ids)
	{ next if (($sample_id eq '.') or ($sample_id eq '..'));
	  my $in_file = "$in_dir/$sample_id/R_SCT/totals.txt";
	  if (!(-e($in_file))) { print "cannot find $in_file\n"; }
	  next if (!(-e($in_file)));
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $NUM_CELLS_KB_UNFILTERED = $line[0]; my $NUM_CELLS_KB_FILTERED = $line[1]; my $NUM_CELLS_WITH_LE_1000_GENES = $line[2]; my $NUM_CELLS_WITH_GE_10k_GENES = $line[3]; my $NUM_CELLS_WITH_LE_2k_UMIS = $line[4]; my $NUM_CELLS_WITH_GE_50k_UMIS = $line[5]; my $NUM_CELLS_WITH_PCT_MT_GE_5 = $line[6]; my $NUM_CELLS_PASSING_ALL_QC_FILTERS = $line[7]; my $NUM_GENES_KB_FILTERED = $line[8]; my $NUM_GENES_PASSING_ALL_QC_FILTERS_AND_IN_GE_3_CELLS = $line[9];
		  my $pc_cells_with_le_1000_genes     = sprintf("%.2f",($NUM_CELLS_WITH_LE_1000_GENES/$NUM_CELLS_KB_FILTERED)*100);
		  my $pc_cells_with_ge_10000_genes    = sprintf("%.2f",($NUM_CELLS_WITH_GE_10k_GENES/$NUM_CELLS_KB_FILTERED)*100);
		  my $pc_cells_with_le_2000_umis      = sprintf("%.2f",($NUM_CELLS_WITH_LE_2k_UMIS/$NUM_CELLS_KB_FILTERED)*100);
		  my $pc_cells_with_ge_50000_umis     = sprintf("%.2f",($NUM_CELLS_WITH_GE_50k_UMIS/$NUM_CELLS_KB_FILTERED)*100);
		  my $pc_cells_with_pct_mt_ge_5       = sprintf("%.2f",($NUM_CELLS_WITH_PCT_MT_GE_5/$NUM_CELLS_KB_FILTERED)*100);
		  my $pc_cells_passing_all_qc_filters = sprintf("%.2f",($NUM_CELLS_PASSING_ALL_QC_FILTERS/$NUM_CELLS_KB_FILTERED)*100);
		  
		  # for each sample, how many cells and genes/cell are in the integrated dataset?
		  my $num_cells_in_integrated_dataset = scalar keys %{$genes_per_sample{$sample_id}};
		  if ($num_cells_in_integrated_dataset != $NUM_CELLS_PASSING_ALL_QC_FILTERS)
			{ print "ERROR: the per-sample Seurat output for $sample_id says that $NUM_CELLS_PASSING_ALL_QC_FILTERS cells pass QC, but the atlas metadata file tells me that this sample contributes $num_cells_in_integrated_dataset\n"; exit 1; }
		  my @num_genes = ();
		  while((my $barcode,my $irrel)=each(%{$genes_per_sample{$sample_id}}))
			{ my $num_genes = $genes_per_sample{$sample_id}{$barcode};
			  push(@num_genes,$num_genes);
			}
		  my @sorted_num_genes = sort {$a <=> $b} @num_genes;
		  my $min_genes_per_cell_in_integrated_dataset	  = $sorted_num_genes[0];
		  my $max_genes_per_cell_in_integrated_dataset	  = $sorted_num_genes[$#sorted_num_genes];
		  my $mean_genes_per_cell_in_integrated_dataset   = avg(@num_genes);    $mean_genes_per_cell_in_integrated_dataset   = ceil($mean_genes_per_cell_in_integrated_dataset);
		  my $median_genes_per_cell_in_integrated_dataset = median(@num_genes); $median_genes_per_cell_in_integrated_dataset = ceil($median_genes_per_cell_in_integrated_dataset);
	  
		  print OUT "$sample_id\t$publications{$sample_id}\t$NUM_CELLS_KB_UNFILTERED\t$NUM_CELLS_KB_FILTERED\t$NUM_CELLS_WITH_LE_1000_GENES\t$pc_cells_with_le_1000_genes\t$NUM_CELLS_WITH_GE_10k_GENES\t$pc_cells_with_ge_10000_genes\t$NUM_CELLS_WITH_LE_2k_UMIS\t$pc_cells_with_le_2000_umis\t$NUM_CELLS_WITH_GE_50k_UMIS\t$pc_cells_with_ge_50000_umis\t$NUM_CELLS_WITH_PCT_MT_GE_5\t$pc_cells_with_pct_mt_ge_5\t$NUM_CELLS_PASSING_ALL_QC_FILTERS\t$pc_cells_passing_all_qc_filters\t$NUM_GENES_PASSING_ALL_QC_FILTERS_AND_IN_GE_3_CELLS\t$min_genes_per_cell_in_integrated_dataset\t$mean_genes_per_cell_in_integrated_dataset\t$median_genes_per_cell_in_integrated_dataset\t$max_genes_per_cell_in_integrated_dataset\n";
		}
	  close(IN) or die $!;
	}
close(OUT) or die $!;
exit 1;