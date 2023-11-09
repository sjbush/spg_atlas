=head

PURPOSE: what proportion of germ cells at each time point are classified as state 0, state 0A, state 0B, and so on?

AFTER USAGE:

library(tidyverse)
theme_set(theme_bw())

df<-read.table('C:/Users/clme2012/Desktop/human_adult_SSCs/states0to4.onto_which_projected_mouse_time_course.proporp_of_cells_mapping_to_each_cluster.min_conf0.8.txt',sep='\t',header=T)
df$Age <- factor(df$Age, levels = c("5-6 days", "10 days", "14-18 days", "20-35 days", "56-67 days", "3-4 months", "5+ months"))
df$Cell.cluster <- factor(df$Cell.cluster, levels = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene"))
#df.sub<-subset(df,df$Total.no..of.cells.with.projection.confidence.score...0.6 > 1000) # & ((df$Cell.cluster == 'state 0') | (df$Cell.cluster == 'state 0A/1') | (df$Cell.cluster == 'state 0B')))
#df<-df.sub

ggplot(df, aes(fill=Cell.cluster, x=as.factor(Age), y=X..of.cells.projected.to.cluster)) + geom_bar(position="fill", stat="identity") + xlab('Age (days)') + ylab('Proportion of cells projected to cluster') + labs(fill='Cell cluster')

ggplot(df, aes(fill=Cell.cluster, x=as.factor(Age), y=No..of.cells.projected.to.cluster)) + geom_bar(position="dodge", stat="identity") + xlab('Age') + ylab('Total no. of cells projected to cluster') + labs(fill='Cell cluster')

=cut

use strict;
use warnings;

# REQUIREMENTS
my $in_file1  = 'states0to4.onto_which_projected_mouse.ERS3000379.metadata.txt';
my $in_file2  = 'states0to4.onto_which_projected_mouse.ERS3000380.metadata.txt';
my $in_file3  = 'states0to4.onto_which_projected_mouse.SRS3990943.metadata.txt';
my $in_file4  = 'states0to4.onto_which_projected_mouse.ERS2575682.metadata.txt';
my $in_file5  = 'states0to4.onto_which_projected_mouse.SRS3990942.metadata.txt';
my $in_file6  = 'states0to4.onto_which_projected_mouse.ERS2575686.metadata.txt';
my $in_file7  = 'states0to4.onto_which_projected_mouse.SRS3990944.metadata.txt';
my $in_file8  = 'states0to4.onto_which_projected_mouse.SRS3990945.metadata.txt';
my $in_file9  = 'states0to4.onto_which_projected_mouse.ERS2575683.metadata.txt';
my $in_file10 = 'states0to4.onto_which_projected_mouse.SRS3990946.metadata.txt';
my $in_file11 = 'states0to4.onto_which_projected_mouse.ERS2575687.metadata.txt';
my $in_file12 = 'states0to4.onto_which_projected_mouse.SRS3990947.metadata.txt';
my $in_file13 = 'states0to4.onto_which_projected_mouse.ERS2575684.metadata.txt';
my $in_file14 = 'states0to4.onto_which_projected_mouse.ERS2575685.metadata.txt';
my $in_file15 = 'states0to4.onto_which_projected_mouse.SRS3990948.metadata.txt';
my $in_file16 = 'states0to4.onto_which_projected_mouse.SRS3990949.metadata.txt';
my $in_file17 = 'states0to4.onto_which_projected_mouse.ERS2575688.metadata.txt';
my $in_file18 = 'states0to4.onto_which_projected_mouse.ERS2575689.metadata.txt';
my $in_file19 = 'states0to4.onto_which_projected_mouse.ERS2575690.metadata.txt';
my $in_file20 = 'states0to4.onto_which_projected_mouse.ERS2575678.metadata.txt';
my $in_file21 = 'states0to4.onto_which_projected_mouse.ERS2575679.metadata.txt';
my $in_file22 = 'states0to4.onto_which_projected_mouse.SRS3097514.metadata.txt';
my $in_file23 = 'states0to4.onto_which_projected_mouse.SRS3097515.metadata.txt';
my $in_file24 = 'states0to4.onto_which_projected_mouse.ERS2575676.metadata.txt';
my $in_file25 = 'states0to4.onto_which_projected_mouse.ERS2575677.metadata.txt';
my $in_file26 = 'states0to4.onto_which_projected_mouse.ERS2575680.metadata.txt';
my $in_file27 = 'states0to4.onto_which_projected_mouse.ERS2575681.metadata.txt';
my $in_file28 = 'states0to4.onto_which_projected_mouse.SRS3189008.metadata.txt';
my $in_file29 = 'states0to4.onto_which_projected_mouse.SRS3189011.metadata.txt';
my $in_file30 = 'states0to4.onto_which_projected_mouse.SRS3189007.metadata.txt';
my $in_file31 = 'states0to4.onto_which_projected_mouse.SRS3189006.metadata.txt';
my $in_file32 = 'states0to4.onto_which_projected_mouse.SRS3189009.metadata.txt';
my $in_file33 = 'states0to4.onto_which_projected_mouse.SRS3189026.metadata.txt';
my $in_file34 = 'states0to4.onto_which_projected_mouse.SRS3189010.metadata.txt';
my $in_file35 = 'states0to4.onto_which_projected_mouse.SRS3189005.metadata.txt';
my $in_file36 = 'states0to4.onto_which_projected_mouse.SRS3189013.metadata.txt';
my $in_file37 = 'states0to4.onto_which_projected_mouse.SRS3189012.metadata.txt';
my @in_files  = ($in_file1,$in_file2,$in_file3,$in_file4,$in_file5,$in_file6,$in_file7,$in_file8,$in_file9,$in_file10,$in_file11,$in_file12,$in_file13,$in_file14,$in_file15,$in_file16,$in_file17,$in_file18,$in_file19,$in_file20,$in_file21,$in_file22,$in_file23,$in_file24,$in_file25,$in_file26,$in_file27,$in_file28,$in_file29,$in_file30,$in_file31,$in_file32,$in_file33,$in_file34,$in_file35,$in_file36,$in_file37);
my $fatal     = 0;
foreach my $in_file (@in_files)
	{ if (!(-e($in_file)))
		{ print "ERROR: unable to find $in_file\n";
		  $fatal++;
		}
	}
exit 1 if ($fatal > 0);

# PARAMETERS
my $min_conf = 0.6;
my %ages = (
'SRS3097510' => '56-67 days',
'SRS3097511' => '56-67 days',
'SRS3097512' => '56-67 days',
'SRS3097513' => '56-67 days',
'SRS3097514' => '56-67 days',
'SRS3097515' => '56-67 days',
'ERS2575676' => '56-67 days',
'ERS2575677' => '56-67 days',
'ERS2575678' => '56-67 days',
'ERS2575679' => '56-67 days',
'ERS2575680' => '56-67 days',
'ERS2575681' => '56-67 days',
'ERS2575682' => '10 days',
'ERS2575683' => '20-35 days',
'ERS2575684' => '20-35 days',
'ERS2575685' => '20-35 days',
'ERS2575686' => '14-18 days',
'ERS2575687' => '20-35 days',
'ERS2575688' => '56-67 days',
'ERS2575689' => '56-67 days',
'ERS2575690' => '56-67 days',
'ERS3000379' => '5-6 days',
'ERS3000380' => '5-6 days',
'SRS3990943' => '5-6 days',
'SRS3990942' => '14-18 days',
'SRS3990944' => '14-18 days',
'SRS3990945' => '14-18 days',
'SRS3990946' => '20-35 days',
'SRS3990947' => '20-35 days',
'SRS3990948' => '56-67 days',
'SRS3990949' => '56-67 days',
'SRS3189004' => '3-4 months',
'SRS3189005' => '5+ months',
'SRS3189006' => '3-4 months',
'SRS3189007' => '3-4 months',
'SRS3189008' => '3-4 months',
'SRS3189009' => '3-4 months',
'SRS3189010' => '5+ months',
'SRS3189011' => '3-4 months',
'SRS3189012' => '5+ months',
'SRS3189013' => '5+ months',
'SRS3189026' => '3-4 months'
);
my %days = (
'SRS3097510' => 63,
'SRS3097511' => 63,
'SRS3097512' => 63,
'SRS3097513' => 63,
'SRS3097514' => 63,
'SRS3097515' => 63,
'ERS2575676' => 67,
'ERS2575677' => 67,
'ERS2575678' => 63,
'ERS2575679' => 63,
'ERS2575680' => 67,
'ERS2575681' => 67,
'ERS2575682' => 10,
'ERS2575683' => 20,
'ERS2575684' => 30,
'ERS2575685' => 35,
'ERS2575686' => 15,
'ERS2575687' => 25,
'ERS2575688' => 61,
'ERS2575689' => 61,
'ERS2575690' => 62,
'ERS3000379' => 5,
'ERS3000380' => 5,
'SRS3990943' => 6,
'SRS3990942' => 14,
'SRS3990944' => 18,
'SRS3990945' => 18,
'SRS3990946' => 25,
'SRS3990947' => 30,
'SRS3990948' => 56,
'SRS3990949' => 56,
'SRS3189004' => 111,
'SRS3189005' => 157,
'SRS3189006' => 120,
'SRS3189007' => 117,
'SRS3189008' => 91,
'SRS3189009' => 120,
'SRS3189010' => 155,
'SRS3189011' => 92,
'SRS3189012' => 284,
'SRS3189013' => 231,
'SRS3189026' => 122
);

# OUTPUT
my $out_file1 = "states0to4.onto_which_projected_mouse_time_course.proporp_of_cells_mapping_to_each_cluster.full_table.min_conf$min_conf.txt";
my $out_file2 = "states0to4.onto_which_projected_mouse_time_course.proporp_of_cells_mapping_to_each_cluster.min_conf$min_conf.txt";
open(OUT1,'>',$out_file1) or die $!; open(OUT2,'>',$out_file2) or die $!;
print OUT1 "Sample ID\tAge\tAge (days)\tTotal no. of cells with projection confidence score > $min_conf\t";
print OUT1 "No. of cells STATE 0\t% of cells STATE 0\t";
print OUT1 "No. of cells STATE 0A/1\t% of cells STATE 0A/1\t";
print OUT1 "No. of cells STATE 0B\t% of cells STATE 0B\t";
print OUT1 "No. of cells SPG\t% of cells SPG\t";
print OUT1 "No. of cells DIFF SPG\t% of cells DIFF SPG\t";
print OUT1 "No. of cells LEPTOTENE\t% of cells LEPTOTENE\t";
print OUT1 "No. of cells ZYGOTENE\t% of cells ZYGOTENE\n";
print OUT2 "Sample ID\tAge\tTotal no. of cells with projection confidence score > $min_conf\tNo. of cells projected to cluster\t% of cells projected to cluster\tCell cluster\n";

my %cells_per_cluster = ();
foreach my $in_file (@in_files)
	{ my $sample_id;
	  if ($in_file =~ /^.+mouse\.(.*?)\.metadata\.txt$/)
		{ $sample_id = $1; }
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  my $predicted_cluster_score = $line[23]; my $predicted_cluster = $line[24];
		  next if ($predicted_cluster_score < $min_conf);
		  $cells_per_cluster{$sample_id}{$predicted_cluster}++;
		}
	  close(IN) or die $!;
	}
my @sample_ids = ();
while((my $sample_id,my $irrel)=each(%cells_per_cluster))
	{ push(@sample_ids,$sample_id); }
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
foreach my $sample_id (@sorted_sample_ids)
	{ if (!(exists($ages{$sample_id}))) { print "ERROR: no age found for $sample_id\n"; exit 1; }
	  my $age  = $ages{$sample_id};
	  my $days = $days{$sample_id};
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
	  print OUT1 "$sample_id\t$age\t$days\t$total\t";
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