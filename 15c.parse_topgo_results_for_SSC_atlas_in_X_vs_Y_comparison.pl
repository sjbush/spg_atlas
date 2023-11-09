use strict;
use warnings;
use LWP::Simple;

# PARAMETERS
my $resolution = 1.2; # 1.1; # 0.3;
my $min_logfc  = 0;
my $total_no_of_genes_with_this_go_term = 0;
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
my $no_of_comparisons = 0;
foreach my $group1 (@group1)
	{ foreach my $group2 (@group2)
		{ next if ($group1 eq $group2);
		  $no_of_comparisons++;
		}
	}

# OUTPUT
my $out_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.go_enrichment_per_cluster.resolution_$resolution".".$no_of_comparisons"."_X_vs_Y_comparisons".".onlyPosTRUE_minPCT25_logFC25.min_log_fc_$min_logfc.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "Comparison: genes differentially (and positively) expressed in state X relative to state Y\tGO category\tGO ID\tGO term description\tTotal no. of genes annotated with this GO term\tObserved no. of genes with this GO term\tExpected no. of genes with this GO term\tFold increase between obs. and exp.\t% of the total genes with this GO term present in the cluster\tp\n";

foreach my $group1 (@group1)
	{ foreach my $group2 (@group2)
		{ next if ($group1 eq $group2);
		  print "$group1, $group2...\n";
		  
		  # REQUIREMENTS
		  my $in_dir = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.go_enrichment_per_cluster.resolution_$resolution".".$group1"."_vs_$group2".".onlyPosTRUE_minPCT25_logFC25.min_log_fc_$min_logfc"; # from 14c.run_topgo_on_cluster_markers_in_X_vs_Y_comparison.pl
		  if (!(-d($in_dir))) { print "ERROR: cannot find $in_dir\n"; exit 1; }

		  opendir(DIR,$in_dir) or die $!;
		  my @subdirs = readdir(DIR);
		  closedir(DIR) or die $!;
		  my @cluster_ids = ();
		  foreach my $cluster_id (@subdirs)
			{ next if (($cluster_id eq '.') or ($cluster_id eq '..'));
			  for(my $x=0;$x<=2;$x++)
				{ my $in_file; my $go_category;
				  if 	($x == 0) { $in_file = "$in_dir/$cluster_id/bp.txt"; $go_category = 'biological process'; }
				  elsif ($x == 1) { $in_file = "$in_dir/$cluster_id/mf.txt"; $go_category = 'molecular function'; }
				  elsif ($x == 2) { $in_file = "$in_dir/$cluster_id/cc.txt"; $go_category = 'cellular component'; }
				  next if (!(-e($in_file)));
				  next if ($go_category ne 'biological process');
				  open(GO,$in_file) or die $!;
				  while(<GO>)
					{ next if ($. == 1);
					  my $go_line = $_; chomp($go_line);
					  my @go_line = split(/\t/,$go_line);
					  my $go_id   = $go_line[0];
					  my $go_term = $go_line[1];
					  my $total_no_with_this_term = $go_line[2];
					  my $obs_with_this_term 	  = $go_line[3];
					  my $exp_with_this_term 	  = $go_line[4];
					  my $go_p 					  = $go_line[5];
					  my $pc_with_this_term 	  = sprintf("%.2f",(($obs_with_this_term/$total_no_with_this_term)*100));
					  next if ($exp_with_this_term == 0);
					  my $fold_increase = sprintf("%.2f",($obs_with_this_term/$exp_with_this_term));
#					  next if ($fold_increase < 2); # FILTER: if the observed number of genes with this GO term in the cluster exceeds the expected number by fewer than 2-fold
					  next if ($total_no_with_this_term < $total_no_of_genes_with_this_go_term); # FILTER: GO terms represented by too few a set of genes
					  next if (($go_p =~ /^\d+/) && ($go_p > 0.05)); # FILTER: GO terms insignificantly represented in this cluster				  
					  print OUT "$group1 relative to $group2\t$go_category\t$go_id\t$go_term\t$total_no_with_this_term\t$obs_with_this_term\t$exp_with_this_term\t$fold_increase\t$pc_with_this_term\t$go_p\n";
					}
				  close(GO) or die $!;
				}
			}
		}
	}
close(OUT) or die $!;
exit 1;