use strict;
use warnings;

# REQUIREMENTS
my $species = 'Homo_sapiens';
my $go_cats = "/project/GorielyLab2021/sbush/ssc_atlas/go_terms.$species.txt"; # from BioMart (Ens106/GRCh38.p13): Gene stable ID, Gene name, GO term evidence code, GO term accession
if (!(-e($go_cats))) { print "ERROR: cannot find $go_cats\n"; exit 1; }

# PARAMETERS
my $min_no_of_genes_per_tissue = 10;

# OUTPUT
my $mappings1 = "/project/GorielyLab2021/sbush/ssc_atlas/go_categories.by_ens_id.$species.map";
my $mappings2 = "/project/GorielyLab2021/sbush/ssc_atlas/go_categories.by_name.$species.map";
open(OUT_MAP1,'>',$mappings1) or die $!; open(OUT_MAP2,'>',$mappings2) or die $!;

# CONVERT GO TERM LIST INTO A FORMAT SUITABLE FOR THE R PACKAGE TOPGO
my %go_cats_by_id = (); my %go_cats_by_name = ();
open(IN,$go_cats) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0]; my $gene_name = $line[1];
	  next if (!(defined($line[2])));
	  my $evidence_code = $line[2]; my $go_term_accession = $line[3];
	  next if (($evidence_code eq 'NAS') or ($evidence_code eq 'ND')); # non-traceable author statement, or no biological data available
	  $go_cats_by_id{$gene_id}{$go_term_accession}++;
	  $go_cats_by_name{$gene_name}{$go_term_accession}++ unless ($gene_name eq '');	  
	}
close(IN) or die $!;
for(my $x=0;$x<=1;$x++)
	{ my %go_cats = ();
	  if 	($x == 0) { %go_cats = %go_cats_by_id;   }
	  elsif ($x == 1) { %go_cats = %go_cats_by_name; }
	  while((my $gene_id,my $irrel)=each(%go_cats))
		{ my $go_term_line = '';
		  while((my $go_term,my $irrel)=each(%{$go_cats{$gene_id}}))
			{ $go_term_line .= "$go_term, "; }
		  $go_term_line =~ s/\, $//;
		  if    ($x == 0) { print OUT_MAP1 "$gene_id\t$go_term_line\n"; }
		  elsif ($x == 1) { print OUT_MAP2 "$gene_id\t$go_term_line\n"; }
		}
	}
close(OUT_MAP1) or die $!; close(OUT_MAP2) or die $!;

exit 1;