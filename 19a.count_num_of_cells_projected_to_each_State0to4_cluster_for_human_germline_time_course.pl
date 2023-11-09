=head

PURPOSE: what proportion of germ cells at each time point are classified as state 0, state 0A, state 0B, and so on?

AFTER USAGE:

library(tidyverse)
theme_set(theme_bw())

df<-read.table('C:/Users/clme2012/Desktop/human_adult_SSCs/states0to4.onto_which_projected_human_germline_time_course.proporp_of_cells_mapping_to_each_cluster.min_conf0.8.txt',sep='\t',header=T)
#df$Age <- factor(df$Age, levels = c("E 6 weeks", "E 8 weeks", "E 9 weeks", "E 10 weeks", "E 15 weeks", "E 19 weeks", "E 23 weeks", "2 days", "5 months", "7 days", "1 year", "13 months", "2 years", "5 years", "7 years", "8 years", "11 years", "13 years", "14 years", "17 years", "20-25 years", "22 years", "24 years", "25 years", "30-40 years", "31 years", "33 years", "34 years", "36 years", "37 years", "42 years", "49 years", "50 years", "55 years", "62 years", "66 years"))
df$Age <- factor(df$Age, levels = c("foetal (8-23 weeks)", "0-2 years", "5-8 years", "11 years", "13-25 years", "30-40 years", "> 40 years"))
df$Cell.cluster <- factor(df$Cell.cluster, levels = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene"))
#df.sub<-subset(df,df$Total.no..of.cells.with.projection.confidence.score...0.8 > 1000) # & ((df$Cell.cluster == 'state 0') | (df$Cell.cluster == 'state 0A/1') | (df$Cell.cluster == 'state 0B')))
#df<-df.sub

ggplot(df, aes(fill=Cell.cluster, x=as.factor(Age), y=X..of.cells.projected.to.cluster)) + geom_bar(position="fill", stat="identity") + xlab('Age') + ylab('Proportion of cells projected to cluster') + labs(fill='Cell cluster')

ggplot(df, aes(fill=Cell.cluster, x=as.factor(Age), y=No..of.cells.projected.to.cluster)) + geom_bar(position="dodge", stat="identity") + xlab('Age') + ylab('Total no. of cells projected to cluster') + labs(fill='Cell cluster')

=cut

use strict;
use warnings;

# PARAMETERS
my $min_conf = 0;
my %ages = (
'HRR131888' => 'foetal (8-23 weeks)',
'HRR131892' => 'foetal (8-23 weeks)',
'HRR131894' => 'foetal (8-23 weeks)',
'HRR131895' => 'foetal (8-23 weeks)',
'HRR131897' => 'foetal (8-23 weeks)',
'HRR131901' => 'foetal (8-23 weeks)',
'HRR131904' => 'foetal (8-23 weeks)',
'HRR131906' => 'foetal (8-23 weeks)',
'SRS4181123' => '0-2 years',
'SRS4181126' => '0-2 years',
'SRS7727466' => '0-2 years',
'SRS7727467' => '0-2 years',
'SRS12015709' => '0-2 years',
'SRS3822680' => '0-2 years',
'SRS3822682' => '0-2 years',
'SRS3822683' => '0-2 years',
'SRS3822686' => '0-2 years',
'LZ011' => '0-2 years',
'SRS12015710' => '0-2 years',
'LZ009' => '5-8 years',
'SRS12015708' => '5-8 years',
'SRS5086057' => '5-8 years',
'SRS5086058' => '5-8 years',
'LZ005' => '5-8 years',
'LZ008' => '11-25 years',
'SRS5086059' => '11-25 years',
'SRS5086060' => '11-25 years',
'SRS5086061' => '11-25 years',
'SRS5086062' => '11-25 years',
'SRS5086063' => '11-25 years',
'SRS5086064' => '11-25 years',
'SRS3065428' => '11-25 years',
'LZ016' => '11-25 years',
'SRS9921715' => '11-25 years',
'SRS5883824' => '11-25 years',
'SRS5883825' => '11-25 years',
'SRS5883826' => '11-25 years',
'SRS9921716' => '11-25 years',
'SRS9921717' => '11-25 years',
'SRS3065429' => '11-25 years',
'SRS3065430' => '11-25 years',
'SRS5883812' => '30-40 years',
'SRS5883813' => '30-40 years',
'SRS5883814' => '30-40 years',
'LZ003' => '30-40 years',
'LZ007' => '30-40 years',
'SRS6959440' => '30-40 years',
'SRS6959441' => '30-40 years',
'SRS2823407' => '30-40 years',
'SRS2823409' => '30-40 years',
'SRS5883810' => '30-40 years',
'SRS5883811' => '30-40 years',
'SRS4181127' => '30-40 years',
'SRS4181130' => '> 40 years',
'SRS2823408' => '> 40 years',
'SRS5086065' => '> 40 years',
'SRS5086066' => '> 40 years',
'SRS6959442' => '> 40 years',
'SRS9921722' => '> 40 years',
'SRS9921719' => '> 40 years',
'SRS9921720' => '> 40 years',
'SRS9921721' => '> 40 years',
'SRS9921723' => '> 40 years');

# REQUIRMENTS
my $fatal = 0;
while((my $sample_id,my $irrel)=each(%ages))
	{ my $in_file1 = "adult_whole_testes.onto_which_projected.$sample_id.cells_which_project_onto_germline.txt"; # from 28.project_human_time_course_onto_adult_whole_testes.R
	  if (!(-e($in_file1)))
		{ print "ERROR: unable to find $in_file1\n";
		  $fatal++;
		}
	  my $in_file2 = "states0to4.onto_which_projected_human_germline.$sample_id.metadata.txt"; # from 29.project_human_germline_time_course_onto_States0to4_resolution1_1.R
	  if (!(-e($in_file2)))
		{ print "ERROR: unable to find $in_file2\n";
		  $fatal++;
		}
	}
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file1 = "states0to4.onto_which_projected_human_germline_time_course.proporp_of_cells_mapping_to_each_cluster.full_table.min_conf$min_conf.txt";
my $out_file2 = "states0to4.onto_which_projected_human_germline_time_course.proporp_of_cells_mapping_to_each_cluster.min_conf$min_conf.txt";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT1 "Sample ID\tAge\tTotal no. of cells with projection confidence score > $min_conf\t";
print OUT1 "No. of cells STATE 0\t% of cells STATE 0\t";
print OUT1 "No. of cells STATE 0A/1\t% of cells STATE 0A/1\t";
print OUT1 "No. of cells STATE 0B\t% of cells STATE 0B\t";
print OUT1 "No. of cells SPG\t% of cells SPG\t";
print OUT1 "No. of cells DIFF SPG\t% of cells DIFF SPG\t";
print OUT1 "No. of cells LEPTOTENE\t% of cells LEPTOTENE\t";
print OUT1 "No. of cells ZYGOTENE\t% of cells ZYGOTENE\n";
print OUT2 "Sample ID\tAge\tTotal no. of cells with projection confidence score > $min_conf\tNo. of cells projected to cluster\t% of cells projected to cluster\tCell cluster\n";

my %cells_per_cluster = ();
while((my $sample_id,my $irrel)=each(%ages))
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
	  my $in_file2 = "states0to4.onto_which_projected_human_germline.$sample_id.metadata.txt"; # from 29.project_human_germline_time_course_onto_States0to4_resolution1_1.R
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
while((my $sample_id,my $irrel)=each(%ages))
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
	{ if (!(exists($ages{$sample_id}))) { print "ERROR: no age found for $sample_id\n"; exit 1; }
	  my $age  		  = $ages{$sample_id};
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
	  print OUT1 "$sample_id\t$age\t$total\t";
	  print OUT1 "$num_state0\t$pc_state0\t";
	  print OUT1 "$num_state0a\t$pc_state0a\t";
	  print OUT1 "$num_state0b\t$pc_state0b\t";
	  print OUT1 "$num_spg\t$pc_spg\t";
	  print OUT1 "$num_diff_spg\t$pc_diff_spg\t";
	  print OUT1 "$num_lepto\t$pc_lepto\t";
	  print OUT1 "$num_zygo\t$pc_zygo\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_state0\t$pc_state0\tstate 0\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_state0a\t$pc_state0a\tstate 0A/1\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_state0b\t$pc_state0b\tstate 0B\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_spg\t$pc_spg\tSPG\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_diff_spg\t$pc_diff_spg\tdiff SPG\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_lepto\t$pc_lepto\tleptotene\n";
	  print OUT2 "$sample_id\t$age\t$total\t$num_zygo\t$pc_zygo\tzygotene\n";
	}

close(OUT1) or die $!; close(OUT2) or die $!;
exit 1;