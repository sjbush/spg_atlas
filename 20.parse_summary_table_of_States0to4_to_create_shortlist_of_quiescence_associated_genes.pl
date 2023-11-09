# parse the atlas table to create a shortlist of quiescence-associated genes
# many of these genes are involved in the positive regulation of quiescence in HSCs and are taken from Li 2011, "Quiescence regulators for hematopoietic stem cell": https://www.sciencedirect.com/science/article/pii/S0301472X11000178
# positive regulators of quiescence are those which decrease quiescence when knocked out - cross-reference the Li 2011 table with Table 4 of Nakamuru-Ishizu 2014: https://journals.biologists.com/dev/article/141/24/4656/46542/The-analysis-roles-and-regulation-of-quiescence-in
# note some name changes relative to the Li 2011 table: TP53 (p53), TAL1 (Scl), CDKN1A (p21), CDKN1C (p57), FBXW7 (FBX7), MYC (c-Myc), NR4A2 (Nurr1), SH2B3 (Lnk), STAT5A (STAT5), RB1 (RB), STK11 (LKB1), CHD4 (Mi-2b)

use strict;
use warnings;


# REQUIREMENTS
my $in_file1 = '/project/GorielyLab2021/sbush/human_adult_SSCs/regulators_of_quiescence.txt'; # manually created by reference to Li 2011 and Nakamuru-Ishizu 2014 (see above)
my $in_file2 = '/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.atlas_resolution1.1.txt'; # from 12c.create_summary_table_of_States0to4_resolution_1.1.pl
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }

# OUTPUT
my $out_file = '/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.quiescence_associated_genes_in_atlas_resolution1.1.txt'; 
open(OUT,'>',$out_file) or die $!;
print OUT "Gene\tSource\tType of regulator\t% of state 0 SSCs in which this gene is detected\t% of state 0A/1 SSCs in which this gene is detected\t% of state 0B SSCs in which this gene is detected\tCluster(s) in which this gene is differentially expressed\n";

# STORE THE NAMES OF INTRINSIC REGULATORS OF QUIESCENCE
my %regulators = ();
open(IN,$in_file1) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r//g;
	  my @line = split(/\t/,$line);
	  my $gene = $line[0]; my $source = $line[1]; my $type = $line[2];
	  $regulators{$gene}{source} = $source;
	  $regulators{$gene}{type}   = $type;
	}
close(IN) or die $!;

my %seen = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $pct_state0 = $line[30]; my $pct_state0a = $line[31]; my $pct_state0b = $line[32]; my $clusters_in_which_de = $line[45];
	  if (exists($regulators{$gene_name}))
		{ print OUT "$gene_name\t$regulators{$gene_name}{source}\t$regulators{$gene_name}{type}\t$pct_state0\t$pct_state0a\t$pct_state0b\t$clusters_in_which_de\n";
		  $seen{$gene_name}++;
		}
	}
close(IN) or die $!;
while((my $gene_name,my $irrel)=each(%regulators))
	{ if (!(exists($seen{$gene_name})))
		{ print "did not see $gene_name\n"; }
	}

close(OUT) or die $!;
exit 1;