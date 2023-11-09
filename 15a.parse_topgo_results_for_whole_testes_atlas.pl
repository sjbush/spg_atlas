use strict;
use warnings;
use LWP::Simple;

# PARAMETERS
my $min_logfc = 0; # 1;
my $total_no_of_genes_with_this_go_term = 0;

# REQUIREMENTS
my $in_dir = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.go_enrichment_per_cluster.min_log_fc_$min_logfc"; # from 14a.run_topgo_on_cluster_markers_for_whole_testes_atlas.pl
if (!(-d($in_dir))) { print "ERROR: cannot find $in_dir\n"; exit 1; }

# OUTPUT
my $out_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.go_enrichment_per_cluster.min_log_fc_$min_logfc.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "Cluster\tGO category\tGO ID\tGO term description\tTotal no. of genes annotated with this GO term\tObserved no. of genes with this GO term\tExpected no. of genes with this GO term\tFold increase between obs. and exp.\t% of the total genes with this GO term present in the cluster\tp\n"; # alternative column 3: "GO ID (depth in GO tree, i.e. maximum distance from a parent term)"

opendir(DIR,$in_dir) or die $!;
my @subdirs = readdir(DIR);
closedir(DIR) or die $!;
my @cluster_ids = ();
foreach my $cluster_id (@subdirs)
	{ next if (($cluster_id eq '.') or ($cluster_id eq '..'));
	  push(@cluster_ids,$cluster_id);		  
	}
my @sorted_cluster_ids = sort {$a <=> $b} @cluster_ids;
my $cluster_ids_seen = 0; my $cluster_ids_total = @sorted_cluster_ids;
foreach my $cluster_id (@sorted_cluster_ids)
	{ $cluster_ids_seen++;
	  my $pc = sprintf("%.2f",(($cluster_ids_seen/$cluster_ids_total)*100)); print "GO: $pc%\n";
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
#			  next if ($fold_increase < 2); # FILTER: if the observed number of genes with this GO term in the cluster exceeds the expected number by fewer than 2-fold
			  next if ($total_no_with_this_term < $total_no_of_genes_with_this_go_term); # FILTER: GO terms represented by too few a set of genes
			  next if (($go_p =~ /^\d+/) && ($go_p > 0.05)); # FILTER: GO terms insignificantly represented in this cluster				  
			  
=cut
			  # what is the 'depth' of the GO ID, i.e. furthest distance from a parent, in the GO tree?
			  my $go_digits;
			  if ($go_id =~ /^GO\:(\d+)$/) { $go_digits = $1; }
			  my $url = "http://supfam.cs.bris.ac.uk/SUPERFAMILY/cgi-bin/go.cgi?search=GO%3A"."$go_digits";
			  my $content = get($url);
			  my $data_available = 0;
			  my @distances_to_parent = ();
			  if (defined($content))
				{ my @content = split(/\n/,$content);
				  foreach my $line (@content)
					{ if ($line =~ /\<tr\>\<td align\=right\>\<font color\=\#FF0000\>/)
						{ my @result = split(/\:/,$line);
						  my $partition = $result[0];
						  if ($partition =~ /^.*?(\d+)$/)
							{ my $distance_to_parent = $1;
							  push(@distances_to_parent,$distance_to_parent);
							  $data_available++;
							}
						}
					}
				}
			  my $furthest_distance_to_parent = '.';
			  if ($data_available > 0)
				{ my @sorted_distances_to_parent = sort {$b <=> $a} @distances_to_parent;
				  $furthest_distance_to_parent = $sorted_distances_to_parent[0];
				}
			  my $go_id_with_parent_depth = "$go_id ($furthest_distance_to_parent)";
=cut
			  print OUT "$cluster_id\t$go_category\t$go_id\t$go_term\t$total_no_with_this_term\t$obs_with_this_term\t$exp_with_this_term\t$fold_increase\t$pc_with_this_term\t$go_p\n";
			}
		  close(GO) or die $!;
		}
	}
close(OUT) or die $!;
exit 1;