=head

BEFORE USAGE:

module add R-base/4.1.0
module add R-cbrg/202111 # R-base/4.1.0 gsl/2.6 hdf5/1.10.7

=cut

use strict;
use warnings;

# PARAMETERS
my $min_cluster_size = 0;
my $min_logfc = 0; # 1;

# REQUIREMENTS
my $clusters = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.cluster_markers.txt"; # from 6.refine_clustering_of_whole_testes_atlas.R
my $mappings = "/project/GorielyLab2021/sbush/human_adult_SSCs/go_categories.by_name.Homo_sapiens.map"; # from 13.turn_go_terms_into_map_file.pl
my $fatal    = 0;
if (!(-e($clusters))) { $fatal++; print "ERROR: cannot find $clusters\n"; }
if (!(-e($mappings))) { $fatal++; print "ERROR: cannot find $mappings\n"; }
if ($fatal > 0) { exit 1; }

# OUTPUT
my $out_dir 			  = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.go_terms_per_cluster.min_log_fc_$min_logfc";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
my $out_dir_for_go_script = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.go_enrichment_per_cluster.min_log_fc_$min_logfc";
my $r_script 			  = "/project/GorielyLab2021/sbush/human_adult_SSCs/run_topgo_analysis_min_log_fc_$min_logfc.R";
open(OUT_R,'>',$r_script) or die $!;
print OUT_R "library(topGO)\n";

# STORE, PER GENE, THE CLUSTER IT BELONGS TO. WE DO THIS ONLY FOR SIGNIFICANTLY DIFFERENTIALLY EXPRESSED GENES WITH A POSITIVE LOG FOLD CHANGE, I.E. GENES THAT ARE *EXPLICITLY CHARACTERISTIC* OF THIS CLUSTER *COMPARED TO ALL OTHERS*
my %clusters = ();
open(IN,$clusters) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/[\r\n]//;
	  my @line = split(/\t/,$line);
	  my $avg_logfc = $line[1]; my $pct_1 = $line[2]; my $pct_2 = $line[3]; my $adj_p = $line[4]; my $cluster_id = $line[5]; my $gene_id = $line[6];
	  $cluster_id =~ s/\//\-/; $cluster_id =~ s/ /\_/; # to prevent downstream errors when it comes to creating $out_dir/$cluster_id.tsv
	  if (($avg_logfc > $min_logfc) and ($adj_p < 0.05))
		{ $clusters{$cluster_id}{$gene_id}++; }
#	  my $abs_diff = abs($pct_1-$pct_2);
#	  if (($avg_logfc > 1) and ($pct_1 > 0.25) and ($adj_p < 0.05))
#		{ $clusters{$cluster_id}{$gene_id}++; }
	}
close(IN) or die $!;
my @cluster_ids = ();
while((my $cluster_id,my $irrel)=each(%clusters))
	{ push(@cluster_ids,$cluster_id); }
my @sorted_cluster_ids = sort {$a cmp $b } @cluster_ids;
foreach my $cluster_id (@sorted_cluster_ids)
	{ my $size = scalar keys %{$clusters{$cluster_id}};
	  next if ($size < $min_cluster_size);
	  open(OUT,'>',"$out_dir/$cluster_id.tsv") or die $!;
	  print OUT "GeneID\n";
	  while((my $gene_id,my $irrel)=each(%{$clusters{$cluster_id}}))
		{ print OUT "$gene_id\n"; }
	  close(OUT) or die $!;
	}

# CREATE R SCRIPT OF TOPGO CODE SO AS TO RUN A GO TERM ENRICHMENT ANALYSIS ON EACH COEXPRESSION CLUSTER
opendir(DIR,$out_dir) or die $!;
my @cluster_files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_cluster_files = sort {$a cmp $b} @cluster_files;
foreach my $cluster_file (@sorted_cluster_files)
	{ next if (($cluster_file eq '.') or ($cluster_file eq '..'));
	  my $cluster_dir;
	  if ($cluster_file =~ /^(.*?)\.tsv$/) { $cluster_dir = $1; }
	  
	  # create empty output directories to be filled later, by this script
	  if (!(-d("$out_dir_for_go_script"))) 			    { mkdir "$out_dir_for_go_script" 			  or die $!; }
	  if (!(-d("$out_dir_for_go_script/$cluster_dir"))) { mkdir "$out_dir_for_go_script/$cluster_dir" or die $!; }
	  print OUT_R "table<-read.table('$out_dir/$cluster_file',sep='\\t',header=T)\n";
	  print OUT_R "cluster_genes<-table\$GeneID\n";
	  print OUT_R "geneID2GO <- readMappings(file = file('$mappings'))\n";
	  print OUT_R "geneNames <- names(geneID2GO)\n";
	  print OUT_R "geneList <- factor(as.integer(geneNames %in% cluster_genes))\n";
	  print OUT_R "names(geneList) <- geneNames\n";
#	  if (!(-e("$out_dir_for_go_script/$cluster_dir/mf.txt")))
#		{ print OUT_R "GOdata <- new('topGOdata', ontology = 'MF', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
#		  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
#		  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
#		  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
#		  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/$cluster_dir/mf.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
#		}
	  if (!(-e("$out_dir_for_go_script/$cluster_dir/bp.txt")))
		{ print OUT_R "GOdata <- new('topGOdata', ontology = 'BP', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
		  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
		  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
		  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
		  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/$cluster_dir/bp.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
		}
#	  if (!(-e("$out_dir_for_go_script/$cluster_dir/cc.txt")))
#		{ print OUT_R "GOdata <- new('topGOdata', ontology = 'CC', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
#		  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
#		  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
#		  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
#		  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/$cluster_dir/cc.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
#		}
	}
close(OUT_R) or die $!;
exit 1;