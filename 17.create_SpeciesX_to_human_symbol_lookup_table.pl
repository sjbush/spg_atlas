use strict;
use warnings;

# REQUIREMENTS
my $species = 'mouse'; # 'rat'; # 'yak'; # 'baboon'; # 'chicken'; # 'panda'; # 'rhesus_macaque'; # 'cynomolgus_macaque'; # 'pig'; # 'sheep';
my $in_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Ens108.$species"."_orthologues_of_human_genes.txt"; # from BioMart: Gene stable ID, Gene name, X gene stable ID, X gene name, X homology type, X orthology confidence [0 low, 1 high
if (!(-e($in_file))) { print "ERROR: cannot find $in_file\n"; exit 1; }

# OUTPUT
my $out_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/$species"."_to_human_symbol_lookup.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "Symbol\tHumanSymbol\n";

my %one_to_ones = ();
open(IN,$in_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[2]))); # CHECKPOINT: there is no point continuing but there is no X gene associated with this human gene
	  my $human_gene_id = $line[0]; my $human_gene_name = $line[1];
	  if ($human_gene_name eq '') { $human_gene_name = $human_gene_id; }
	  my $X_gene_id = $line[2]; my $X_gene_name = $line[3];
	  if ($X_gene_name eq '') { $X_gene_name = $X_gene_id; }
	  my $orthology_type = $line[4];
	  if ($orthology_type eq 'ortholog_one2one')
		{ $one_to_ones{$X_gene_name}{$human_gene_name}++; }
	}
close(IN) or die $!;
my @gene_names = ();
while((my $X_gene_name,my $irrel)=each(%one_to_ones))
	{ push(@gene_names,$X_gene_name); }
my @sorted_gene_names = sort {"$a" cmp "$b"} @gene_names;
foreach my $X_gene_name (@sorted_gene_names)
	{ my $num_orthologues = scalar keys %{$one_to_ones{$X_gene_name}};
	  next if ($num_orthologues != 1);
	  while((my $human_gene_name,my $irrel)=each(%{$one_to_ones{$X_gene_name}}))
		{ print OUT "$X_gene_name\t$human_gene_name\n";
		}
	}

close(OUT) or die $!;
exit 1;