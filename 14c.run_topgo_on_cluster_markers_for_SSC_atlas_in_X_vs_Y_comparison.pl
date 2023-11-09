=head

BEFORE USAGE:

module add R-base/4.1.0
module add R-cbrg/202111 # R-base/4.1.0 gsl/2.6 hdf5/1.10.7

=cut

use strict;
use warnings;

# PARAMETERS
my $min_cluster_size = 0;
my $min_logfc 		 = 0; # 1;
my $resolution 		 = 1.2; # 1.1; # 0.3;
my @group1 = (); my @group2 = ();
if ($resolution == 0.3)
	{ @group1 = (qw/Adark Apale/);
	  @group2 = (qw/Adark Apale/);
	}
elsif ($resolution == 1.1)
	{ @group1 = (qw/State0 State0A State0B SPG meiosis/);
	  @group2 = (qw/State0 State0A State0B SPG meiosis/);
	}
elsif ($resolution == 1.2)
	{ @group1 = (qw/State0 State0A State0to0B State0Bto1/);
  	  @group2 = (qw/State0 State0A State0to0B State0Bto1/);
	}

# OUTPUT	
my $r_script = "/project/GorielyLab2021/sbush/human_adult_SSCs/run_topgo_analyses_resolution_$resolution.onlyPosTRUE_minPCT25_logFC25.min_log_fc_$min_logfc.R";
open(OUT_R,'>',$r_script) or die $!;
print OUT_R "library(topGO)\n";
	
foreach my $group1 (@group1)
	{ foreach my $group2 (@group2)
		{ next if ($group1 eq $group2);
		  print "$group1, $group2...\n";

		  # REQUIREMENTS
		  my $clusters = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.degs_"."$group1"."_vs_"."$group2"."_resolution_$resolution"."_onlyPosTRUE_minPCT25_logFC25.txt"; # from 10.run_DEG_for_States0to4.R
		  my $mappings = "/project/GorielyLab2021/sbush/human_adult_SSCs/go_categories.by_name.Homo_sapiens.map"; # from 13.turn_go_terms_into_map_file.pl
		  my $fatal    = 0;
		  if (!(-e($clusters))) { $fatal++; print "ERROR: cannot find $clusters\n"; }
		  if (!(-e($mappings))) { $fatal++; print "ERROR: cannot find $mappings\n"; }
		  if ($fatal > 0) { exit 1; }

		  # OUTPUT
		  my $out_dir 			    = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.go_terms_per_cluster.resolution_$resolution".".$group1"."_vs_$group2".".onlyPosTRUE_minPCT25_logFC25.min_log_fc_$min_logfc";
		  my $out_dir_for_go_script = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.go_enrichment_per_cluster.resolution_$resolution".".$group1"."_vs_$group2".".onlyPosTRUE_minPCT25_logFC25.min_log_fc_$min_logfc";
		  if (!(-d($out_dir))) { mkdir $out_dir or die $!; }
		  
		  # STORE THE SET OF DE GENES FOR THIS SINGULAR, CLUSTER 1 VS. CLUSTER 2, COMPARISON
		  my %clusters = ();
		  open(IN,$clusters) or die $!;
		  while(<IN>)
			{ next if ($. == 1);
			  my $line = $_; chomp($line); $line =~ s/[\r\n]//;
			  my @line = split(/\t/,$line);
			  my $gene_id = $line[0]; my $avg_logfc = $line[2]; my $adj_p = $line[5];
			  if (($avg_logfc > $min_logfc) and ($adj_p < 0.05))
				{ $clusters{$group1}{$gene_id}++; }
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
			  if (!(-e("$out_dir_for_go_script/$cluster_dir/bp.txt")))
				{ print OUT_R "GOdata <- new('topGOdata', ontology = 'BP', allGenes = geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)\n";
				  print OUT_R "test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test')\n";
				  print OUT_R "resultWeight <- getSigGroups(GOdata, test.stat)\n";
				  print OUT_R "allRes <- GenTable(GOdata, weight = resultWeight)\n";
				  print OUT_R "write.table(allRes,file='$out_dir_for_go_script/$cluster_dir/bp.txt',quote = FALSE, sep = '\\t', row.names = FALSE, col.names = TRUE)\n";
				}
			}
		}
	}
close(OUT_R) or die $!;
exit 1;