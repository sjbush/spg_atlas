=head

PURPOSE: what proportion of germ cells for each species are classified as state 0, state 0A, state 0B, and so on?

AFTER USAGE:

library(tidyverse)
theme_set(theme_bw())

df<-read.table('C:/Users/clme2012/Desktop/human_adult_SSCs/states0to4.onto_which_projected_other_species.proporp_of_cells_mapping_to_each_cluster.min_conf0.8.txt',sep='\t',header=T)
df$Desc <- factor(df$Age, levels = c("sheep, 1.5 years old", "pig, 150 days old", "buffalo, 3 months old", "buffalo, 2 years old", "cynomolgus macaque, 1 year", "cynomolgus macaque, 2 years", "cynomolgus macaque, > 4 years (Lau 2020 'adult 1')", "cynomolgus macaque, > 4 years (Lau 2020 'adult 2')"))
df$Cell.cluster <- factor(df$Cell.cluster, levels = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene"))

ggplot(df, aes(fill=Cell.cluster, x=as.factor(Age), y=X..of.cells.projected.to.cluster)) + geom_bar(position="fill", stat="identity") + xlab('Age') + ylab('Proportion of cells projected to cluster') + labs(fill='Cell cluster')

ggplot(df, aes(fill=Cell.cluster, x=as.factor(Age), y=No..of.cells.projected.to.cluster)) + geom_bar(position="dodge", stat="identity") + xlab('Age') + ylab('Total no. of cells projected to cluster') + labs(fill='Cell cluster')

=cut

use strict;
use warnings;

# PARAMETERS
my $min_conf = 0.8;
my %samples = (
'SRS7528324' => 'sheep, 1.5 years old',
'SRS9029393' => 'pig, 150 days old',
'SRS11257878' => 'buffalo, 3 months old',
'SRS11257879' => 'buffalo, 2 years old',
'infant' => 'cynomolgus macaque, 1 year',
'juvenile' => 'cynomolgus macaque, 2 years',
'adult1' => "cynomolgus macaque, > 4 years (Lau 2020 'adult 1')",
'adult2' => "cynomolgus macaque, > 4 years (Lau 2020 'adult 2')",
'SRS9839846' => 'rat, 8-10 weeks',
'SRS9839847' => 'rat, 8-10 weeks',
'SRS9839848' => 'rat, 8-10 weeks',
'SRS9839849' => 'rat, 8-10 weeks',
'SRS9839850' => 'rat, 8-10 weeks',
'SRS9839851' => 'rat, 8-10 weeks',
'SRS9839852' => 'rat, 8-10 weeks',
'SRS9839853' => 'rat, 8-10 weeks'
);
my %species = (
'SRS7528324' => 'sheep',
'SRS9029393' => 'pig',
'SRS11257878' => 'buffalo',
'SRS11257879' => 'buffalo',
'infant' => 'cynomolgus_macaque',
'juvenile' => 'cynomolgus_macaque',
'adult1' => 'cynomolgus_macaque',
'adult2' => 'cynomolgus_macaque',
'SRS9839846' => 'rat',
'SRS9839847' => 'rat',
'SRS9839848' => 'rat',
'SRS9839849' => 'rat',
'SRS9839850' => 'rat',
'SRS9839851' => 'rat',
'SRS9839852' => 'rat',
'SRS9839853' => 'rat'
);

# REQUIREMENTS
my $fatal = 0;
while((my $sample_id,my $species)=each(%species))
	{ my $in_file1 = "adult_whole_testes.onto_which_projected.$sample_id.cells_which_project_onto_germline.txt"; # from 30.project_sheep_pig_buffalo_and_cynomolgus_macaque_onto_adult_whole_testes.R and 40.project_rat_onto_adult_whole_testes.R 
	  if (!(-e($in_file1)))
		{ print "ERROR: unable to find $in_file1\n";
		  $fatal++;
		}
	  my $in_file2 = "states0to4.onto_which_projected_$species"."_germline.$sample_id.metadata.txt"; # from 31.project_sheep_pig_buffalo_and_cynomolgus_macaque_onto_States0to4_resolution1_1.R and 41.project_rat_onto_States0to4_resolution1_1.R
	  if (!(-e($in_file2)))
		{ print "ERROR: unable to find $in_file2\n";
		  $fatal++;
		}
	}
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file1 = "states0to4.onto_which_projected_other_species.proporp_of_cells_mapping_to_each_cluster.full_table.min_conf$min_conf.txt";
my $out_file2 = "states0to4.onto_which_projected_other_species.proporp_of_cells_mapping_to_each_cluster.min_conf$min_conf.txt";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT1 "Sample ID\tDescription\tTotal no. of cells with projection confidence score > $min_conf\t";
print OUT1 "No. of cells STATE 0\t% of cells STATE 0\t";
print OUT1 "No. of cells STATE 0A/1\t% of cells STATE 0A/1\t";
print OUT1 "No. of cells STATE 0B\t% of cells STATE 0B\t";
print OUT1 "No. of cells SPG\t% of cells SPG\t";
print OUT1 "No. of cells DIFF SPG\t% of cells DIFF SPG\t";
print OUT1 "No. of cells LEPTOTENE\t% of cells LEPTOTENE\t";
print OUT1 "No. of cells ZYGOTENE\t% of cells ZYGOTENE\n";
print OUT2 "Sample ID\tDescription\tTotal no. of cells with projection confidence score > $min_conf\tNo. of cells projected to cluster\t% of cells projected to cluster\tCell cluster\n";

my %cells_per_cluster = ();
while((my $sample_id,my $species)=each(%species))
	{ # first identify those cells which map to the germline
	  my %germline_cells = ();
	  my $in_file1 = "adult_whole_testes.onto_which_projected.$sample_id.cells_which_project_onto_germline.txt"; # from 28.project_human_time_course_onto_adult_whole_testes.R
	  open(IN,$in_file1) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $barcode = $line[22];
		  $germline_cells{$barcode}++;
		}
	  close(IN) or die $!;
	  
	  # let's confirm these cells are represented by a unique barcode
	  while((my $barcode,my $num)=each(%germline_cells))
		{ if ($num > 1)
			{ print "ERROR in $sample_id: $num cells have the barcode $barcode\n"; exit 1; }
		}
	  
	  # now identify which cluster each cell projects to, restricting analysis only to germline cells
	  my $in_file2 = "states0to4.onto_which_projected_$species"."_germline.$sample_id.metadata.txt"; # from 29.project_human_germline_time_course_onto_States0to4_resolution1_1.R
	  open(IN,$in_file2) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $barcode = $line[22]; my $predicted_cluster_score = $line[23]; my $predicted_cluster = $line[24];
		  next if (!(exists($germline_cells{$barcode}))); # CHECKPOINT: exclude cells which did not map to the germline when projected onto the whole-testes atlas
		  next if ($predicted_cluster_score < $min_conf); # CHECKPOINT: exclude cells which do meet the minimum cluster score threshold ($min_conf)
		  $cells_per_cluster{$sample_id}{$predicted_cluster}++;
		}
	  close(IN) or die $!;
	}
	
# were there any samples where no cells met the minimum cluster score threshold?
while((my $sample_id,my $irrel)=each(%samples))
	{ if (!(exists($cells_per_cluster{$sample_id})))
		{ print "WARNING: no cells could be included from sample $sample_id because none met the minimum score threshold ($min_conf)\n"; }
	}

my @sample_ids = ();
while((my $sample_id,my $irrel)=each(%cells_per_cluster))
	{ push(@sample_ids,$sample_id); }
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
my $num_samples = @sorted_sample_ids;
print "data was obtained from $num_samples samples\n";
foreach my $sample_id (@sorted_sample_ids)
	{ if (!(exists($samples{$sample_id}))) { print "ERROR: no description found for $sample_id\n"; exit 1; }
	  my $desc  	   = $samples{$sample_id};
	  my $num_state0   = 0; if (exists($cells_per_cluster{$sample_id}{'state 0'}))    { $num_state0   = $cells_per_cluster{$sample_id}{'state 0'};    }
	  my $num_state0a  = 0; if (exists($cells_per_cluster{$sample_id}{'state 0A/1'})) { $num_state0a  = $cells_per_cluster{$sample_id}{'state 0A/1'}; }
	  my $num_state0b  = 0; if (exists($cells_per_cluster{$sample_id}{'state 0B'}))   { $num_state0b  = $cells_per_cluster{$sample_id}{'state 0B'};   }
	  my $num_spg      = 0; if (exists($cells_per_cluster{$sample_id}{'SPG'}))    	  { $num_spg      = $cells_per_cluster{$sample_id}{'SPG'};        }
	  my $num_diff_spg = 0; if (exists($cells_per_cluster{$sample_id}{'diff SPG'}))   { $num_diff_spg = $cells_per_cluster{$sample_id}{'diff SPG'};   }
	  my $num_lepto    = 0; if (exists($cells_per_cluster{$sample_id}{'leptotene'}))  { $num_lepto    = $cells_per_cluster{$sample_id}{'leptotene'};  }
	  my $num_zygo     = 0; if (exists($cells_per_cluster{$sample_id}{'zygotene'}))   { $num_zygo     = $cells_per_cluster{$sample_id}{'zygotene'};   }
	  my $total 	   = $num_state0+$num_state0a+$num_state0b+$num_spg+$num_diff_spg+$num_lepto+$num_zygo;
	  my $pc_state0    = sprintf("%.2f",(($num_state0/$total)*100));
	  my $pc_state0a   = sprintf("%.2f",(($num_state0a/$total)*100));
	  my $pc_state0b   = sprintf("%.2f",(($num_state0b/$total)*100));
	  my $pc_spg       = sprintf("%.2f",(($num_spg/$total)*100));
	  my $pc_diff_spg  = sprintf("%.2f",(($num_diff_spg/$total)*100));
	  my $pc_lepto     = sprintf("%.2f",(($num_lepto/$total)*100));
	  my $pc_zygo      = sprintf("%.2f",(($num_zygo/$total)*100));
	  print OUT1 "$sample_id\t$desc\t$total\t";
	  print OUT1 "$num_state0\t$pc_state0\t";
	  print OUT1 "$num_state0a\t$pc_state0a\t";
	  print OUT1 "$num_state0b\t$pc_state0b\t";
	  print OUT1 "$num_spg\t$pc_spg\t";
	  print OUT1 "$num_diff_spg\t$pc_diff_spg\t";
	  print OUT1 "$num_lepto\t$pc_lepto\t";
	  print OUT1 "$num_zygo\t$pc_zygo\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_state0\t$pc_state0\tstate 0\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_state0a\t$pc_state0a\tstate 0A/1\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_state0b\t$pc_state0b\tstate 0B\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_spg\t$pc_spg\tSPG\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_diff_spg\t$pc_diff_spg\tdiff SPG\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_lepto\t$pc_lepto\tleptotene\n";
	  print OUT2 "$sample_id\t$desc\t$total\t$num_zygo\t$pc_zygo\tzygotene\n";
	}

close(OUT1) or die $!; close(OUT2) or die $!;
exit 1;