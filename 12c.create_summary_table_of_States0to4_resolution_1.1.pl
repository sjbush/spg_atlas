=head
PURPOSE: process the various atlas metadata files to create a summary table that details, per gene, its average expression across all cells per state, and the proportion of cells in which it is detectably expressed in a given state. We cross-reference this information with various sources (OMIM, Mendeliome, DDD, COSMIC, Kaplanis 2020, Seurat differential expression analysis) to predict candidate genes for selfish selection.
=cut

use strict;
use warnings;
use Acme::Tools qw(avg);

# PARAMETERS
my $res = 1.1;
my $total_cells_in_atlas = 4447;
my $total_sscs_in_atlas  = 2229; # i.e. the number of cells in clusters 0, 2 and 4. Count can be obtained from states0to4.metadata.txt, produced by 8.subset_whole_testes_atlas_to_obtain_States0to4.R
my $total_non_sscs_in_atlas = $total_cells_in_atlas-$total_sscs_in_atlas;
my @cluster_order = ("Cluster 2","Cluster 0","Cluster 4","Cluster 1","Cluster 6","Cluster 3","Cluster 5"); # this relates only to the order in which the DE results are presented
my @selfishly_selected_genes = (qw/BRAF CBL FGFR2 FGFR3 HRAS KRAS MAP2K1 PTPN11 RAF1 RET SMAD4 SOS1/);
my @raindance_screened_genes = (qw/ABL1 ACVR1 AKT1 AKT2 AKT3 ALK APC ATP1A3 BRAF CBL CSF1R CTNNB1 DNMT3A EGFR ERBB2 EZH2 FBN1 FGFR1 FGFR2 FGFR3 FLT3 GNAS HRAS IDH1 IDH2 JAK2 KAT6B KIF22 KIT KRAS LRP5 MAP2K1 MAP2K2 MECP2 MED12 MET MTOR NF1 NOTCH1 NOTCH2 NOTCH3 NRAS PDGFRA PIK3CA PIK3R2 PTEN PTPN11 RAF1 RALB RET SETBP1 SHOC2 SKI SMAD4 SOS1 TGFBR1 TGFBR2 TP53 TP63 VHL ZBTB16/); # 61 genes previously screened using RainDance PCR. This list is derived from Supplementary Tables 2 and 6 of Maher 2018. Excluded are 10 'neutral test' negative control genes: BFSP1 BFSP2 CHIC1 COX15 GADPH HMBS MDH1 MIP RHO RPL13A
my %known_selfish 	   = map {$_ => 1} @selfishly_selected_genes;
my %already_raindanced = map {$_ => 1} @raindance_screened_genes;
my $cluster0_name = 'STATE 0A/1';
my $cluster1_name = 'SPG';
my $cluster2_name = 'STATE 0';
my $cluster3_name = 'LEPTOTENE';
my $cluster4_name = 'STATE 0B';
my $cluster5_name = 'ZYGOTENE';
my $cluster6_name = 'DIFF SPG';
my %cluster_names = (
'Cluster 0' => 'state 0A/1',
'Cluster 1' => 'SPG',
'Cluster 2' => 'state 0',
'Cluster 3' => 'leptotene',
'Cluster 4' => 'state 0B',
'Cluster 5' => 'zygotene',
'Cluster 6' => 'diff SPG'
);

# REQUIREMENTS
my $pheno_file 	   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Ens104.gene_annotations.txt'; 		   # from Ensembl BioMart: Gene stable ID, Gene name, Gene description, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Strand, Gene type, Phenotype description, Source name
my $id_to_reactome = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Ens105.geneID_to_reactomeID.txt'; 	   # from BioMart: gene stable ID, Reactome ID
my $orthologues    = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Ens107.mouse_orthologues.txt'; 		   # from BioMart: Gene stable ID, Mouse gene stable ID, Mouse gene name, Last common ancestor with Mouse, Mouse homology type, %id. target Mouse gene identical to query gene, %id. query gene identical to target Mouse gene, Mouse orthology confidence [0 low, 1 high]
my $reactome_descs = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/ReactomePathways.txt'; 				   # from https://reactome.org/download/current/ReactomePathways.txt, accessed 14th August 2022
my $markers  	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.degs_All_vs_All_resolution_$res.txt"; # from 10.run_DE_on_States0to4.R
my $genelist1      = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/OMIM_dominant_inheritance.txt'; 	   # 1828 genes; from "genemap2_22_12_21 KW.xlsx"
my $genelist2  	   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/TruSight_Mendeliome_panel.txt'; 	   # 6669 genes; from "genemap2_22_12_21 KW.xlsx"
my $genelist3  	   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/DDG2P_v3.1.tsv';					   # 1922 genes; from https://panelapp.genomicsengland.co.uk/panels/484/ | https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC4392068/
my $genelist4  	   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/cancer_gene_census.txt'; 			   # 736 genes, irrespective of 'hallmark' or 'tier' status; from https://cancer.sanger.ac.uk/census
my $genelist5  	   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Kaplanis285.txt'; 			  		   # 285 genes; from https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2832-5/MediaObjects/41586_2020_2832_MOESM4_ESM.xlsx and filtered to retain only those where 'diagnostic category' is 'discordant', 'consensus' or 'novel', and 'significance' is 'TRUE' (Kaplanis 2020)
my $genelist6  	   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Kaplanis22.txt';				  	   # 15 genes from https://www.nature.com/articles/s41586-020-2832-5/tables/1 + 7 from Spencer
my $t2g  	   	   = '/project/GorielyLab2021/sbush/ssc_atlas/indexes/Homo_sapiens/t2g2.txt';				   		   # from 3.update_t2g_to_replace_null_names.pl
my $proporp1   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster0_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $proporp2   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster1_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $proporp3   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster2_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $proporp4   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster3_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $proporp5   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster4_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $proporp6   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster5_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $proporp7   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.proporp_of_cells_expressing_each_gene_in_Cluster6_resolution_$res.txt";  # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num1   	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster0_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num2   	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster1_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num3  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster2_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num4  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster3_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num5  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster4_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num6  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster5_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $num7  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.number_of_cells_expressing_each_gene_in_Cluster6_resolution_$res.txt";   # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr1  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster0_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr2  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster1_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr3  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster2_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr4  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster3_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr5  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster4_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr6  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster5_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $expr7  	   	   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster6_resolution_$res.txt"; # from 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R
my $twist_no_cov   = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Target_bases_not_covered_by_probes_Oxford_800genes_1X_TE-99535631_hg38.bed'; # emailed by TWIST on 11th April
my $fatal 		   = 0;
my @in_files = ($pheno_file,$id_to_reactome,$reactome_descs,$orthologues,$markers,$genelist1,$genelist2,$genelist3,$genelist4,$genelist5,$genelist6,$t2g,$proporp1,$proporp2,$proporp3,$proporp4,$proporp5,$proporp6,$proporp7,$num1,$num2,$num3,$num4,$num5,$num6,$num7,$expr1,$expr2,$expr3,$expr4,$expr5,$expr6,$expr7,$twist_no_cov);
foreach my $in_file (@in_files)
	{ if (!(-e($in_file)))
		{ print "ERROR: unable to find $in_file\n";
		  $fatal++;
		}
	}
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.atlas_resolution$res.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "Gene name (given as name/Ensembl ID if not unique)\tEnsembl gene ID\tGene description\tGene type\tLocation\t";
print OUT "Mouse gene name\tMouse gene ID\tHomology type\t% identity mouse-to-human\t% identity human-to-mouse\tOrthology confidence\t";
print OUT "Phenotype\tPhenotype source\t";
print OUT "Reactome pathway ID(s)\tReactome pathway description(s)\tIs this gene on the oncogenic MAPK signalling pathway? (Reactome ID R-HSA-6802957)\t";
print OUT "Number of CLUSTER 2 ($cluster2_name) cells in which this gene is detected\t";
print OUT "Number of CLUSTER 0 ($cluster0_name) cells in which this gene is detected\t";
print OUT "Number of CLUSTER 4 ($cluster4_name) cells in which this gene is detected\t";
print OUT "Number of CLUSTER 1 ($cluster1_name) cells in which this gene is detected\t";
print OUT "Number of CLUSTER 6 ($cluster6_name) cells in which this gene is detected\t";
print OUT "Number of CLUSTER 3 ($cluster3_name) cells in which this gene is detected\t";
print OUT "Number of CLUSTER 5 ($cluster5_name) cells in which this gene is detected\t";
print OUT "Total number of cells in which this gene is detected\t% of cells in which this gene is detected\t";
print OUT "Total number of undifferentiated SSCs (i.e. CLUSTER 0 + CLUSTER 2 + CLUSTER 4) in which this gene is detected\t% of undifferentiated SSCs (i.e. CLUSTER 0 + CLUSTER 2 + CLUSTER 4) in which this gene is detected\t";
print OUT "Total number of non-SSC cells in which this gene is detected (i.e. all except CLUSTERS 0, 2 and 4)\t% of non-SSC cells in which this gene is detected (i.e. all except CLUSTERS 0, 2 and 4)\t";
print OUT "Ratio of % of SSC to non-SSCs in which this gene is detected (only calculated if gene is present in > 0 non-SSC cells)\t";
print OUT "% of CLUSTER 2 ($cluster2_name) cells in which this gene is detected\t";
print OUT "% of CLUSTER 0 ($cluster0_name) cells in which this gene is detected\t";
print OUT "% of CLUSTER 4 ($cluster4_name) cells in which this gene is detected\t";
print OUT "% of CLUSTER 1 ($cluster1_name) cells in which this gene is detected\t";
print OUT "% of CLUSTER 6 ($cluster6_name) cells in which this gene is detected\t";
print OUT "% of CLUSTER 3 ($cluster3_name) cells in which this gene is detected\t";
print OUT "% of CLUSTER 5 ($cluster5_name) cells in which this gene is detected\t";
print OUT "Avg. expression across all CLUSTER 2 ($cluster2_name) cells\t";
print OUT "Avg. expression across all CLUSTER 0 ($cluster0_name) cells\t";
print OUT "Avg. expression across all CLUSTER 4 ($cluster4_name) cells\t";
print OUT "Avg. expression across all CLUSTER 1 ($cluster1_name) cells\t";
print OUT "Avg. expression across all CLUSTER 6 ($cluster6_name) cells\t";
print OUT "Avg. expression across all CLUSTER 3 ($cluster3_name) cells\t";
print OUT "Avg. expression across all CLUSTER 5 ($cluster5_name) cells\t";
print OUT "Number of clusters in which this gene is differentially expressed\t";
print OUT "Cluster(s)\tProportion of cells expressing this gene in this cluster\tProportion of cells expressing this gene in all other clusters\tAbsolute difference in proportion of cells expressing this gene in this cluster versus all other clusters\tMaximum absolute difference in proportion of cells expressing this gene in this cluster versus all other clusters\tavg. log2 fold change\tp-value\tadj. p-value\t";
print OUT "Cluster for which this gene is a candidate biomarker (criteria: differentially expressed only in this cluster)\t";
print OUT "Is this gene DE in CLUSTER 2 ($cluster2_name)?\t";
print OUT "Is this gene DE in CLUSTER 0 ($cluster0_name)?\t";
print OUT "Is this gene DE in CLUSTER 4 ($cluster4_name)?\t";
print OUT "Is this gene DE in CLUSTER 1 ($cluster1_name)?\t";
print OUT "Is this gene DE in CLUSTER 6 ($cluster6_name)?\t";
print OUT "Is this gene DE in CLUSTER 3 ($cluster3_name)?\t";
print OUT "Is this gene DE in CLUSTER 5 ($cluster5_name)?\t";
print OUT "Is this gene OMIM dominant?\t";
print OUT "Is this gene in the Mendeliome?\t";
print OUT "Is this gene on the DDD panel with status green ('confirmed')?\t";
print OUT "Is this gene on the DDD panel with status green ('confirmed') or amber ('probable')?\t";
print OUT "Is this gene in COSMIC?\t";
print OUT "Is this gene in Kaplanis 2020?\t";
print OUT "No. of gene lists in which this gene appears\t";
print OUT "Is this gene already known to be selfishly selected?\t";
print OUT "Has this gene already been screened using RainDance (in Maher et al. 2018)?\t";
print OUT "Is some region within this gene NOT covered by a TWIST probe?\n";

# STORE HUMAN-TO-MOUSE ORTHOLOGUES
my %orthologues = ();
open(IN,$orthologues) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $human_gene_id = $line[0];
	  next if (!(defined($line[1])));
	  my $mouse_gene_id = $line[1]; my $mouse_gene_name = $line[2]; my $homology_type = $line[4]; my $pc_ident_mouse_to_human = $line[5]; my $pc_ident_human_to_mouse = $line[6]; my $confidence = $line[7]; 
	  if ($mouse_gene_name eq '') { $mouse_gene_name = $mouse_gene_id; }
	  my $orthologue_line = "$mouse_gene_name\t$mouse_gene_id\t$homology_type\t$pc_ident_mouse_to_human\t$pc_ident_human_to_mouse\t$confidence";
	  $orthologues{$human_gene_id} = $orthologue_line;
	}
close(IN) or die $!;

# STORE GENE ANNOTATIONS
my %gene_info = (); my %gene_name_to_gene_id_lookup = ();
open(IN,$pheno_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_id = $line[0];
	  my $gene_name = ''; my $gene_desc = '.'; my $chr = '.'; my $gene_start = '.'; my $gene_end = '.'; my $strand = '.'; my $gene_type = '.'; my $phenotype = '.'; my $source = '.';
	  if (defined($line[1]))  { $gene_name   = $line[1]; }
	  if (defined($line[2]))  { $gene_desc   = $line[2]; }
	  if (defined($line[3]))  { $chr 	     = $line[3]; }
	  if (defined($line[4]))  { $gene_start  = $line[4]; }
	  if (defined($line[5]))  { $gene_end    = $line[5]; }
	  if (defined($line[6]))  { $strand      = $line[6]; }
	  if (defined($line[7]))  { $gene_type   = $line[7]; }
	  if (defined($line[8]))  { $phenotype   = $line[8]; }
	  if (defined($line[9]))  { $source      = $line[9]; }
	  next if ($chr =~ /CHR\_/); # ignore alternative sequence gene IDs
	  if ($gene_name   eq '') { $gene_name = $gene_id; }
	  my $loc = "$chr:$gene_start-$gene_end:$strand";
	  if ($gene_desc   eq '') { $gene_desc   = '.'; }
	  if ($gene_type   eq '') { $gene_type   = '.'; }
	  if ($loc 		   eq '') { $loc 		 = '.'; }
	  if ($phenotype   eq '') { $phenotype   = '.'; }
	  if ($source      eq '') { $source      = '.'; }
	  $gene_name_to_gene_id_lookup{$gene_name}{$gene_id}++;
	  $gene_info{$gene_id}{desc} = $gene_desc;
	  $gene_info{$gene_id}{type} = $gene_type;
	  $gene_info{$gene_id}{loc}  = $loc;
	  $gene_info{$gene_id}{phenot}{$phenotype}++;
	  $gene_info{$gene_id}{source}{$source}++;
	}
close(IN) or die $!;

# STORE OMIM 'DOMINANT PHENOTYPE' GENES
my %all_genes = (); my %omim_genes = ();
open(IN,$genelist1) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0];
	  $omim_genes{$gene_name}++;
	  $all_genes{$gene_name}++;
	}
close(IN) or die $!;

# STORE MENDELIOME GENES
my %mendeliome = ();
open(IN,$genelist2) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0];
	  $mendeliome{$gene_name}++;
	  $all_genes{$gene_name}++;
	}
close(IN) or die $!;

# STORE DDD DATA, TO WHICH WE ALSO ADD THE LONGLIST OF 285 DEVELOPMENTAL-DISORDER ASSOCIATED GENES FROM KAPLANIS, ET AL.
my %ddd_green = (); my %ddd_green_by_name = ();
my %ddd_green_and_amber = (); my %ddd_green_and_amber_by_name = ();
open(IN,$genelist3) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[2]; my $source = $line[3]; my $gene_id = $line[21];
	  if ($source =~ /Expert Review Green/)
		{ $ddd_green{$gene_id}++;
		  $ddd_green_by_name{$gene_name}++;
		  $all_genes{$gene_name}++;
		}
	  if (($source =~ /Expert Review Green/) or ($source =~ /Expert Review Amber/))
		{ $ddd_green_and_amber{$gene_id}++;
		  $ddd_green_and_amber_by_name{$gene_name}++;
		  $all_genes{$gene_name}++;
		}
	}
close(IN) or die $!;
open(IN,$genelist5) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0];
	  $ddd_green_by_name{$gene_name}++;
	  $ddd_green_and_amber_by_name{$gene_name}++;
	  $all_genes{$gene_name}++;
	}
close(IN) or die $!;

# STORE COSMIC TUMORIGENESIS GENES
my %cancer_genes = ();
open(IN,$genelist4) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0];
	  $cancer_genes{$gene_name}++;
	  $all_genes{$gene_name}++;
	}
close(IN) or die $!;

# STORE A SHORTLIST OF 15+7 DEVELOPMENTAL-DISORDER ASSOCIATED GENES FROM KAPLANIS, ET AL. THESE ARE GENES WHICH HAVE RECURRENT DE NOVO SNVs.
my %kaplanis_genes = ();
open(IN,$genelist6) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0];
	  $kaplanis_genes{$gene_name}++;
	  $all_genes{$gene_name}++;
	}
close(IN) or die $!;

# STORE T2G LOOKUP
my %t2g = ();
open(IN,$t2g) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $transcript_id  = $line[0]; my $gene_id = $line[1];
	  if ($transcript_id =~ /(.+?)\.\d+$/) { $transcript_id = $1; }
	  if ($gene_id		 =~ /(.+?)\.\d+$/) { $gene_id 		= $1; }
	  $t2g{$transcript_id} = $gene_id;
	}
close(IN) or die $!;

# STORE THE % OF CELLS EXPRESSING A PARTICULAR GENE IN EACH SSC STATE
my %pct = ();
for(my $x=0;$x<=6;$x++)
	{ my $in_file = ''; my $state = "Cluster $x";
	  if    ($x == 0) { $in_file = $proporp1; }
	  elsif ($x == 1) { $in_file = $proporp2; }
	  elsif ($x == 2) { $in_file = $proporp3; }
	  elsif ($x == 3) { $in_file = $proporp4; }
	  elsif ($x == 4) { $in_file = $proporp5; }
	  elsif ($x == 5) { $in_file = $proporp6; }
	  elsif ($x == 6) { $in_file = $proporp7; }
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/\r$//;
		  my @line = split(/\t/,$line);
		  my $gene_name = $line[0]; my $proportion = $line[1];
		  my $pct = sprintf("%.3f",($proportion*100));
		  $pct{$gene_name}{$state} = $pct;
		}
	  close(IN) or die $!;
	}
	
# STORE THE NUMBER OF CELLS EXPRESSING A PARTICULAR GENE IN EACH SSC STATE
my %ct = ();
for(my $x=0;$x<=6;$x++)
	{ my $in_file = ''; my $state = "Cluster $x";
	  if    ($x == 0) { $in_file = $num1; }
	  elsif ($x == 1) { $in_file = $num2; }
	  elsif ($x == 2) { $in_file = $num3; }
	  elsif ($x == 3) { $in_file = $num4; }
	  elsif ($x == 4) { $in_file = $num5; }
	  elsif ($x == 5) { $in_file = $num6; }
	  elsif ($x == 6) { $in_file = $num7; }
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/\r$//;
		  my @line = split(/\t/,$line);
		  my $gene_name = $line[0]; my $ct = $line[1];
		  $ct{$gene_name}{$state} = $ct;
		}
	  close(IN) or die $!;
	}
	
# STORE THE AVG. EXPRESSION OF EACH GENE IN EACH SSC STATE
my %avg = ();
for(my $x=0;$x<=6;$x++)
	{ my $in_file = ''; my $state = "Cluster $x";
	  if    ($x == 0) { $in_file = $expr1; }
	  elsif ($x == 1) { $in_file = $expr2; }
	  elsif ($x == 2) { $in_file = $expr3; }
	  elsif ($x == 3) { $in_file = $expr4; }
	  elsif ($x == 4) { $in_file = $expr5; }
	  elsif ($x == 5) { $in_file = $expr6; }
	  elsif ($x == 6) { $in_file = $expr7; }
	  open(IN,$in_file) or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line); $line =~ s/\r$//;
		  my @line = split(/\t/,$line);
		  my $gene_name = $line[0];
		  my @expr = ();
		  for(my $x=1;$x<@line;$x++)
			{ my $expr = $line[$x];
			  push(@expr,$expr);
			}
		  my $avg = avg(@expr); $avg = sprintf("%.5f",$avg);
		  $avg{$gene_name}{$state} = $avg;
		}
	  close(IN) or die $!;
	}
	
# STORE THE RESULTS OF THE SEURAT DIFFERENTIAL EXPRESSION ANALYSIS, RETAINING ONLY PUTATIVE POSITIVE BIOMARKERS (i.e. THOSE WHERE EXPRESSION IN *THIS* CLUSTER IS GREATER THAN ALL OTHER CLUSTERS)
my %de_genes = ();
open(IN,$markers) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $p_val = $line[0]; my $avg_log2FC = $line[1]; my $pct_1 = $line[2]; my $pct_2 = $line[3]; my $p_val_adj = $line[4]; my $cluster = $line[5]; my $gene_name = $line[6];
	  next if ($pct_1 < $pct_2);
	  next if ($p_val_adj > 0.05);
	  $avg_log2FC = sprintf("%.3f",$avg_log2FC);
	  $p_val 	  = sprintf("%.3e",$p_val);
	  $p_val_adj  = sprintf("%.3e",$p_val_adj);
	  $de_genes{$gene_name}{"Cluster $cluster"}{pct_1} 	    = $pct_1;
	  $de_genes{$gene_name}{"Cluster $cluster"}{pct_2} 	    = $pct_2;
	  $de_genes{$gene_name}{"Cluster $cluster"}{diff} 	    = sprintf("%.3f",(abs($pct_1-$pct_2)));
	  $de_genes{$gene_name}{"Cluster $cluster"}{p_val} 	    = $p_val;
	  $de_genes{$gene_name}{"Cluster $cluster"}{p_val_adj}  = $p_val_adj;
	  $de_genes{$gene_name}{"Cluster $cluster"}{avg_log2FC} = $avg_log2FC;
	}
close(IN) or die $!;

# STORE GENE ID TO REACTOME ID LOOKUP
my %gene_to_reactome_id = ();
open(IN,$id_to_reactome) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  next if (!(defined($line[1])));
	  my $gene_id = $line[0]; my $reactome_id = $line[1];
	  $gene_to_reactome_id{$gene_id}{$reactome_id}++;
	}
close(IN) or die $!;

# STORE REACTOME ID DESCRIPTIONS
my %reactome_descs = ();
open(IN,$reactome_descs) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  my $reactome_id = $line[0]; my $reactome_desc = $line[1]; my $species = $line[2];
	  next if ($species ne 'Homo sapiens');
	  next if ($reactome_id !~ /^R\-/);
	  $reactome_descs{$reactome_id} = $reactome_desc;
	}
close(IN) or die $!;

# STORE THE TWIST 'NO COVERAGE' PROBES
my %twist_no_cov = ();
open(IN,$twist_no_cov) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene = $line[3];
	  if ($gene =~ /^(.*?)\_ENS.+$/) { $gene = $1; }
	  $twist_no_cov{$gene}++;
	}
close(IN) or die $!;

# PRINT COUNT AND EXPRESSION DATA PER GENE, AS A SUMMARY OF THE CONTENTS OF THE ATLAS, AND DETERMINE WHICH CANDIDATE GENES ARE FOUND IN WHICH STATE
my @gene_names = ();
while((my $gene_name,my $irrel)=each(%pct))
	{ push(@gene_names,$gene_name); }
my @sorted_gene_names = sort {"\L$a" cmp "\L$b"} @gene_names;
foreach my $gene_name (@sorted_gene_names)
	{ my $gene_id = '.';
	  my $gene_desc = '.'; my $gene_type = '.'; my $loc = '.'; my $phenotype = '.'; my $source = '.';
	  if ($gene_name =~ /^.*?\/(ENSG.*?)$/)
		{ $gene_id   = $1;
		  $gene_desc = $gene_info{$gene_id}{desc};
	      $gene_type = $gene_info{$gene_id}{type};
		  $loc 		 = $gene_info{$gene_id}{loc};
		  my @phenot = (); my @source = ();
		  while((my $phenot,my $irrel)=each(%{$gene_info{$gene_id}{phenot}})) { push(@phenot,$phenot); }
		  while((my $source,my $irrel)=each(%{$gene_info{$gene_id}{source}})) { push(@source,$source); }
		  my @sorted_phenot = sort {$a cmp $b} @phenot;
		  my @sorted_source = sort {$a cmp $b} @source;
		  $phenotype = join(" | ",@sorted_phenot);
		  $source    = join(" | ",@sorted_source);
		}
	  else
		{ my $num_gene_ids = scalar keys %{$gene_name_to_gene_id_lookup{$gene_name}};
		  if ($num_gene_ids != 1) { print "error: unable to assign a gene ID to $gene_name as there are $num_gene_ids candidates\n"; }
		  if ($num_gene_ids == 1)
			{ while((my $this_gene_id,my $irrel)=each(%{$gene_name_to_gene_id_lookup{$gene_name}}))
				{ $gene_id   = $this_gene_id;
				  $gene_desc = $gene_info{$this_gene_id}{desc};
				  $gene_type = $gene_info{$this_gene_id}{type};
				  $loc 		 = $gene_info{$this_gene_id}{loc};
				  my @phenot = (); my @source = ();
				  while((my $phenot,my $irrel)=each(%{$gene_info{$gene_id}{phenot}})) { push(@phenot,$phenot); }
				  while((my $source,my $irrel)=each(%{$gene_info{$gene_id}{source}})) { push(@source,$source); }
				  my @sorted_phenot = sort {$a cmp $b} @phenot;
				  my @sorted_source = sort {$a cmp $b} @source;
				  $phenotype = join(" | ",@sorted_phenot);
				  $source    = join(" | ",@sorted_source);
				}
			}
		}
	  next if ($gene_id eq '.'); # CHECKPOINT: we are unable to identify an Ensembl gene ID
	  
	  # is this a candidate gene for selfish selection, on the basis of OMIM, the Mendeliome, DDD, COSMIC and/or Kaplanis?
	  my $omim_status 		= 'no';
	  my $mendeliome_status = 'no';
	  my $ddd_green_status  = 'no';
	  my $ddd_green_and_amber_status = 'no';
	  my $cancer_status     = 'no';
	  my $kaplanis_status   = 'no';
	  my $no_of_gene_lists  = 0;
	  if ((exists($ddd_green{$gene_id})) or (exists($ddd_green{$gene_name})))
		{ $ddd_green_status = 'yes'; $no_of_gene_lists++; }
	  if ((exists($ddd_green_and_amber{$gene_id})) or (exists($ddd_green_and_amber{$gene_name})))
		{ $ddd_green_and_amber_status = 'yes'; $no_of_gene_lists++; }
	  if  (exists($omim_genes{$gene_name}))     { $omim_status 	     = 'yes'; $no_of_gene_lists++; }
	  if  (exists($mendeliome{$gene_name}))     { $mendeliome_status = 'yes'; $no_of_gene_lists++; }
	  if  (exists($cancer_genes{$gene_name}))   { $cancer_status     = 'yes'; $no_of_gene_lists++; }
	  if  (exists($kaplanis_genes{$gene_name})) { $kaplanis_status   = 'yes'; $no_of_gene_lists++; }
	  
	  # parse the Seurat DE analysis to determine which clusters this gene is differentially expressed in, if any
	  my %clusters_this_gene_is_de_in = ();
	  if (exists($de_genes{$gene_name}))
		{ while((my $cluster,my $irrel)=each(%{$de_genes{$gene_name}}))
			{ $clusters_this_gene_is_de_in{$cluster}++;
			}
		}
	  my $number_of_clusters_in_which_this_gene_is_de = scalar keys %clusters_this_gene_is_de_in;
	  
	  # is this gene a candidate biomarker for a particular cluster?
	  my $is_biomarker_for = '.';
	  if ($number_of_clusters_in_which_this_gene_is_de == 1)
		{ while((my $cluster,my $irrel)=each(%{$de_genes{$gene_name}}))
			{ $is_biomarker_for = $cluster_names{$cluster};
			}
		}
		
	  # which pathways are this gene associated with, and are any of them Reactome ID R-HSA-6802957 ('oncogenic MAPK signalling')?
	  my $reactome_ids = '.'; my $reactome_descs = '.';
	  my $is_this_gene_on_the_oncogenic_mapk_pathway = 'no';
	  if (exists($gene_to_reactome_id{$gene_id}))
		{ my @reactome_ids = ();
		  while((my $reactome_id,my $irrel)=each(%{$gene_to_reactome_id{$gene_id}}))
			{ next if (!(exists($reactome_descs{$reactome_id})));
			  push(@reactome_ids,$reactome_id);
			  if ($reactome_id eq 'R-HSA-6802957')
				{ $is_this_gene_on_the_oncogenic_mapk_pathway = 'yes'; }
			}
		   my @sorted_reactome_ids = sort {$a cmp $b} @reactome_ids;
		   $reactome_ids = join(" | ",@sorted_reactome_ids);
	  	   $reactome_descs = '';
		   foreach my $reactome_id (@sorted_reactome_ids)
			{ $reactome_descs .= "$reactome_descs{$reactome_id} | "; }
		   $reactome_ids =~ s/\| $//; $reactome_descs =~ s/\| $//;
		}
	  
	  # is this gene already known to be selfishly selected?
	  my $known_selfish = 'no';
	  if (exists($known_selfish{$gene_name}))
		{ $known_selfish = 'yes'; }
	  
	  # has this gene already been screened for selfish selection using RainDance PCR, as in Maher 2018?
	  my $already_raindanced = 'no';
	  if (exists($already_raindanced{$gene_name}))
		{ $already_raindanced = 'yes'; }
	  
	  # print raw count of cells in which this gene is detected, expression data, and differential expression results per gene
	  my $tot_cells_in_which_gene_detected    = $ct{$gene_name}{'Cluster 0'}+$ct{$gene_name}{'Cluster 1'}+$ct{$gene_name}{'Cluster 2'}+$ct{$gene_name}{'Cluster 3'}+$ct{$gene_name}{'Cluster 4'}+$ct{$gene_name}{'Cluster 5'}+$ct{$gene_name}{'Cluster 6'};
	  my $pct_cells_in_which_gene_detected    = sprintf("%.3f",(($tot_cells_in_which_gene_detected/$total_cells_in_atlas)*100));
	  my $tot_sscs_in_which_gene_detected     = $ct{$gene_name}{'Cluster 0'}+$ct{$gene_name}{'Cluster 2'}+$ct{$gene_name}{'Cluster 4'};
	  my $pct_sscs_in_which_gene_detected     = sprintf("%.3f",(($tot_sscs_in_which_gene_detected/$total_sscs_in_atlas)*100));
	  my $tot_non_sscs_in_which_gene_detected = $ct{$gene_name}{'Cluster 1'}+$ct{$gene_name}{'Cluster 3'}+$ct{$gene_name}{'Cluster 5'}+$ct{$gene_name}{'Cluster 6'};
	  my $pct_non_sscs_in_which_gene_detected = sprintf("%.3f",(($tot_non_sscs_in_which_gene_detected/$total_non_sscs_in_atlas)*100));
	  my $ratio_of_pct_sscs_to_pct_non_sscs   = 'NA';
	  if ($pct_non_sscs_in_which_gene_detected > 0) { $ratio_of_pct_sscs_to_pct_non_sscs = sprintf("%.3f",($pct_sscs_in_which_gene_detected/$pct_non_sscs_in_which_gene_detected)); }
	  my $orthologue_line = ".\t.\t.\t.\t.\t.";
	  my $de_in_cluster0 = 'no'; my $de_in_cluster1 = 'no'; my $de_in_cluster2 = 'no'; my $de_in_cluster3 = 'no'; my $de_in_cluster4 = 'no'; my $de_in_cluster5 = 'no'; my $de_in_cluster6 = 'no';
	  if (exists($orthologues{$gene_id})) 	 { $orthologue_line = $orthologues{$gene_id}; }
	  my $not_covered = 'not assessed';
	  if (exists($twist_no_cov{$gene_name})) { $not_covered = 'yes'; }
	  if ($number_of_clusters_in_which_this_gene_is_de == 0)
		{ print OUT "$gene_name\t$gene_id\t$gene_desc\t$gene_type\t$loc\t$orthologue_line\t$phenotype\t$source\t$reactome_ids\t$reactome_descs\t$is_this_gene_on_the_oncogenic_mapk_pathway\t";
		  print OUT "$ct{$gene_name}{'Cluster 2'}\t$ct{$gene_name}{'Cluster 0'}\t$ct{$gene_name}{'Cluster 4'}\t$ct{$gene_name}{'Cluster 1'}\t$ct{$gene_name}{'Cluster 6'}\t$ct{$gene_name}{'Cluster 3'}\t$ct{$gene_name}{'Cluster 5'}\t";
		  print OUT "$tot_cells_in_which_gene_detected\t$pct_cells_in_which_gene_detected\t";
		  print OUT "$tot_sscs_in_which_gene_detected\t$pct_sscs_in_which_gene_detected\t";
		  print OUT "$tot_non_sscs_in_which_gene_detected\t$pct_non_sscs_in_which_gene_detected\t";
		  print OUT "$ratio_of_pct_sscs_to_pct_non_sscs\t";
		  print OUT "$pct{$gene_name}{'Cluster 2'}\t$pct{$gene_name}{'Cluster 0'}\t$pct{$gene_name}{'Cluster 4'}\t$pct{$gene_name}{'Cluster 1'}\t$pct{$gene_name}{'Cluster 6'}\t$pct{$gene_name}{'Cluster 3'}\t$pct{$gene_name}{'Cluster 5'}\t";
		  print OUT "$avg{$gene_name}{'Cluster 2'}\t$avg{$gene_name}{'Cluster 0'}\t$avg{$gene_name}{'Cluster 4'}\t$avg{$gene_name}{'Cluster 1'}\t$avg{$gene_name}{'Cluster 6'}\t$avg{$gene_name}{'Cluster 3'}\t$avg{$gene_name}{'Cluster 5'}\t";
		  print OUT "$number_of_clusters_in_which_this_gene_is_de\t.\t.\t.\t.\t.\t.\t.\t.\t.\t";
		  print OUT "$de_in_cluster2\t$de_in_cluster0\t$de_in_cluster4\t$de_in_cluster1\t$de_in_cluster6\t$de_in_cluster3\t$de_in_cluster5\t";
		  print OUT "$omim_status\t$mendeliome_status\t$ddd_green_status\t$ddd_green_and_amber_status\t$cancer_status\t$kaplanis_status\t$no_of_gene_lists\t";
		  print OUT "$known_selfish\t$already_raindanced\t$not_covered\n";
		}
	  elsif ($number_of_clusters_in_which_this_gene_is_de > 0)
		{ my $de_line_state = ''; my $de_line_pct_1 = ''; my $de_line_pct_2 = ''; my $de_line_diff = ''; my $de_line_logfc = ''; my $de_line_p_val = ''; my $de_line_adj_p = '';
		  my @diffs = ();
		  foreach my $cluster (@cluster_order)
			{ if (exists($de_genes{$gene_name}{$cluster}))
				{ $de_line_state .= "$cluster ($cluster_names{$cluster}) | ";
				  $de_line_pct_1 .= "$de_genes{$gene_name}{$cluster}{pct_1} | ";
				  $de_line_pct_2 .= "$de_genes{$gene_name}{$cluster}{pct_2} | ";
				  $de_line_diff  .= "$de_genes{$gene_name}{$cluster}{diff} | ";
				  $de_line_logfc .= "$de_genes{$gene_name}{$cluster}{avg_log2FC} | ";
				  $de_line_p_val .= "$de_genes{$gene_name}{$cluster}{p_val} | ";
				  $de_line_adj_p .= "$de_genes{$gene_name}{$cluster}{p_val_adj} | ";
				  push(@diffs,$de_genes{$gene_name}{$cluster}{diff});
				  if 	($cluster eq 'Cluster 0') { $de_in_cluster0 = 'yes'; }
				  elsif ($cluster eq 'Cluster 1') { $de_in_cluster1 = 'yes'; }
				  elsif ($cluster eq 'Cluster 2') { $de_in_cluster2 = 'yes'; }
				  elsif ($cluster eq 'Cluster 3') { $de_in_cluster3 = 'yes'; }
				  elsif ($cluster eq 'Cluster 4') { $de_in_cluster4 = 'yes'; }
				  elsif ($cluster eq 'Cluster 5') { $de_in_cluster5 = 'yes'; }
				  elsif ($cluster eq 'Cluster 6') { $de_in_cluster6 = 'yes'; }
				}
			}
		  $de_line_state =~ s/ \| $//; $de_line_pct_1 =~ s/ \| $//; $de_line_pct_2 =~ s/ \| $//; $de_line_diff =~ s/ \| $//; $de_line_logfc =~ s/ \| $//; $de_line_p_val =~ s/ \| $//; $de_line_adj_p =~ s/ \| $//;
		  my @sorted_diffs = sort {$b <=> $a} @diffs;
		  my $max_diff = $sorted_diffs[0];
		  print OUT "$gene_name\t$gene_id\t$gene_desc\t$gene_type\t$loc\t$orthologue_line\t$phenotype\t$source\t$reactome_ids\t$reactome_descs\t$is_this_gene_on_the_oncogenic_mapk_pathway\t";
		  print OUT "$ct{$gene_name}{'Cluster 2'}\t$ct{$gene_name}{'Cluster 0'}\t$ct{$gene_name}{'Cluster 4'}\t$ct{$gene_name}{'Cluster 1'}\t$ct{$gene_name}{'Cluster 6'}\t$ct{$gene_name}{'Cluster 3'}\t$ct{$gene_name}{'Cluster 5'}\t";
		  print OUT "$tot_cells_in_which_gene_detected\t$pct_cells_in_which_gene_detected\t";
		  print OUT "$tot_sscs_in_which_gene_detected\t$pct_sscs_in_which_gene_detected\t";
		  print OUT "$tot_non_sscs_in_which_gene_detected\t$pct_non_sscs_in_which_gene_detected\t";
		  print OUT "$ratio_of_pct_sscs_to_pct_non_sscs\t";
		  print OUT "$pct{$gene_name}{'Cluster 2'}\t$pct{$gene_name}{'Cluster 0'}\t$pct{$gene_name}{'Cluster 4'}\t$pct{$gene_name}{'Cluster 1'}\t$pct{$gene_name}{'Cluster 6'}\t$pct{$gene_name}{'Cluster 3'}\t$pct{$gene_name}{'Cluster 5'}\t";
		  print OUT "$avg{$gene_name}{'Cluster 2'}\t$avg{$gene_name}{'Cluster 0'}\t$avg{$gene_name}{'Cluster 4'}\t$avg{$gene_name}{'Cluster 1'}\t$avg{$gene_name}{'Cluster 6'}\t$avg{$gene_name}{'Cluster 3'}\t$avg{$gene_name}{'Cluster 5'}\t";
		  print OUT "$number_of_clusters_in_which_this_gene_is_de\t$de_line_state\t$de_line_pct_1\t$de_line_pct_2\t$de_line_diff\t$max_diff\t$de_line_logfc\t$de_line_p_val\t$de_line_adj_p\t$is_biomarker_for\t";
		  print OUT "$de_in_cluster2\t$de_in_cluster0\t$de_in_cluster4\t$de_in_cluster1\t$de_in_cluster6\t$de_in_cluster3\t$de_in_cluster5\t";
		  print OUT "$omim_status\t$mendeliome_status\t$ddd_green_status\t$ddd_green_and_amber_status\t$cancer_status\t$kaplanis_status\t$no_of_gene_lists\t";
		  print OUT "$known_selfish\t$already_raindanced\t$not_covered\n";
		}
	}
close(OUT) or die $!;
exit 1;