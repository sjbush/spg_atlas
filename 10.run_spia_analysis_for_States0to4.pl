=head
BEFORE USAGE:
module add R-cbrg/202111

AFTER USAGE:

SPIA requires a multi-part dataset so that the command spia() can function. Necessities are (1) a vector of the log2 fold changes of the differentially expressed genes, the names of the vector being Entrez gene IDs , (2) a list of the Entrez gene IDs in the reference set, (3) a KEGG organism ID.
However...
Be aware that this warning can occur: “The KEGG pathway data for your organism is not present in the extdata folder of the SPIA package!!! Please generate one first using makeSPIAdata and specify its location using data.dir argument or copy it in the extdata folder of the SPIA package!” (See the manual for details regarding this function: http://www.bioconductor.org/packages/release/bioc/manuals/SPIA/man/SPIA.pdf). This is particularly obvious when attempting to run spia() on non-humans.
The problem is simple: we need to first obtain the necessary data, in the form of KGML (‘KEGG XML’) files (an exchange format for KEGG pathway maps), all which are available from the KEGG website (http://www.kegg.jp/kegg/xml/). You can either download KGML files individually (i.e. for each pathway map) or in bulk by using the ‘get’ method of the KEGG API. Unfortunately, this is paywalled.

Note also that the local install of SPIA/ comes with an extdata/hsa dir pre-loaded with only 4 maps (hsa03013, hsa03050, hsa04914, hsa05210) – as a consequence, the spia() function will work nicely with the example code that comes with the package, although this is human-only and with limited scope. This is useless for data of broader scope and/or encompassing more species. As such, these KGML files need both supplementing and updating or the package is effectively just a useless diversion – but how to do this if you don’t want to pay?

0. First check to see if the species has KEGG data available. See the organism IDs here: http://www.genome.jp/kegg/catalog/org_list.html
1. Go to the list of pathway maps: http://www.genome.jp/kegg/pathway.html. We want to get KGML files corresponding to all of these. Don’t bother doing it manually – just note that each pathway map has a numerical ID, e.g. 00052 = galactose metabolism.
2. Go to the KEGG BRITE database: http://www.genome.jp/kegg/brite.html. This contains links to ‘htext’ files, tab-delimited descriptions of the binary relationships between KEGG objects (such as what we’re after, the pathway maps). Select ‘KEGG pathway maps’ under the ‘Brite mapping: Pathways and Ontologies’ category: http://www.genome.jp/kegg-bin/get_htext?br08901.keg
3. Now select ‘download htext.’ The file you obtain, br08901.keg, now lists the name of each pathway with its associated number.
	--> this file is accurate as of 8th March 2022
4. Each number in this file can be prefixed with a three-letter organism code. If you want to know about, e.g., human galatose metabolism, look up hsa00052 (and you can download the corresponding KGML file accordingly). We’ll use this htext file as the input to a shell script to download all the KGML files we like.
5. Run the following, which will download every pathway.
	--> curl -s http://rest.kegg.jp/list/pathway/hsa | awk '{split($1,a,":"); print "curl -s http://rest.kegg.jp/get/"a[2]"/kgml -o "a[2]".xml"}' | bash # see https://www.r-bloggers.com/2018/06/download-all-kegg-pathway-kgml-files-for-spia-analysis/
6. Note that not all of the output files will contain data; this is because KEGG is not a database complete in all particulars and so certain pathways will not be mapped in particular species. Move all non-empty files into SPIA/extdata/keggxml, subdir ‘hsa’.
7. Finally, from R, run the following: “makeSPIAdata(kgml.path=‘C:/Users/sbush/Documents/R/win-library/4.1/SPIA/extdata/keggxml/hsa’,organism=‘hsa’,out.path=‘C:/Users/sbush/Documents/R/win-library/4.1/SPIA/extdata')”
8. Now to run the SPIA analysis, do this:

library(SPIA)

table<-read.table('C:/Users/sbush/Desktop/human_adult_SSCs/gene_lists/Ens105.ens_to_entrez_ids.txt',sep='\t',header=T)
all_vector<-as.character(table$NCBI.gene..formerly.Entrezgene..ID)
all_vector <- all_vector[!is.na(all_vector)]
all_vector <- unique(all_vector)

table<-read.table('C:/Users/sbush/Desktop/human_adult_SSCs/results/states0to4.degs_State0A_vs_State0_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.forSPIA.txt',sep='\t',header=T)
df<-data.frame(ENTREZ_ID=table$EntrezID,LOG2FC=table$log2FC)
de_vector<-df$LOG2FC
names(de_vector)<-as.character(df$ENTREZ_ID)
res_df<-spia(de=de_vector,all=all_vector,organism='hsa',nB=2000,plots=FALSE,verbose=TRUE,beta=NULL,combine='fisher')
number_of_rows<-nrow(res_df)
out_df <- data.frame()
name<-res_df$Name
id<-res_df$ID
psize<-res_df$pSize
nde<-res_df$NDE
pnde<-res_df$pNDE
ta<-res_df$tA
ppert<-res_df$pPERT
pg<-res_df$pG
pgfdr<-res_df$pGFdr
pgfwer<-res_df$pGFWER
status<-res_df$Status
newrow<-data.frame(NAME=name,ID=id,pSIZE=psize,NDE=nde,pNDE=pnde,tA=ta,pPERT=ppert,pG=pg,pGFdr=pgfdr,STATUS=status)
out_df=rbind(out_df,newrow)
write.table(out_df,'C:/Users/sbush/Desktop/human_adult_SSCs/results/states0to4.degs_State0A_vs_State0_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.kegg_enrichment_for_state_0_DEGs.txt',quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

9. visualise log fold changes on a pathway as so:

# the headers for this output file are as follows: KEGG pathway name\tKEGG pathway ID\tTotal number of genes on this pathway\tNumber of differentially expressed genes on this pathway\tProbability of observing at least this number of differentially expressed genes on this pathway by chance\tTotal accumulated perturbation of this pathway (TA)\tProbability of observing a total perturbation of this pathway greater than TA by chance\tGlobal probability value\tGlobal probability value (FDR 5%)\tPathway status

library(dplyr)
library(pathview)

df<-read.table('C:/Users/sbush/Desktop/human_adult_SSCs/Ens105_ens_to_entrez_ids.txt',sep='\t',header=T)
name_to_entrez<-data.frame(NAME=df$Gene.name,ENT=df$NCBI.gene..formerly.Entrezgene..ID)

df<-read.table('C:/Users/sbush/Desktop/human_adult_SSCs/SSC_STATES_0_TO_2_contents_of_SSC_atlas.txt',sep='\t',header=T,quote="")
#logfc<-data.frame(NAME=df$Gene.name..given.as.name.Ensembl.ID.if.not.unique.,LOG_FC=log2(df$X..of.STATE.0B.cells.in.which.this.gene.is.detected/df$X..of.STATE.0A.cells.in.which.this.gene.is.detected))
logfc<-data.frame(NAME=df$Gene.name..given.as.name.Ensembl.ID.if.not.unique.,LOG_FC=log2(df$Avg..expression.across.all.STATE.0B.cells/df$Avg..expression.across.all.STATE.0A.cells))
logfc <- do.call(data.frame, lapply(logfc, function(x) replace(x, is.infinite(x), NA))) # see https://statisticsglobe.com/replace-inf-with-na-in-vector-and-data-frame-in-r
logfc<-na.omit(logfc)
merged_df<-merge(name_to_entrez,logfc,by='NAME')
merged_df<-na.omit(merged_df)
exp<-merged_df$LOG_FC
gene_id<-merged_df$ENT
exp_per_gene<-data.frame(EXP=exp,ID=gene_id)
exp_per_gene<-exp_per_gene[!duplicated(exp_per_gene$ID),]
exp_per_gene<-data.frame(exp_per_gene$EXP,row.names=exp_per_gene$ID)
exp_per_gene<-na.omit(exp_per_gene)

# look up pathway IDs here: https://www.genome.jp/kegg/pathway.html
pathview(gene.data=exp_per_gene, pathway.id = "05012", species = "hsa", out.suffix = "parkinson")
pathview(gene.data=exp_per_gene, pathway.id = "05016", species = "hsa", out.suffix = "huntington")
pathview(gene.data=exp_per_gene, pathway.id = "04110", species = "hsa", out.suffix = "cell_cycle")
pathview(gene.data=exp_per_gene, pathway.id = "05168", species = "hsa", out.suffix = "herpes_simplex")
pathview(gene.data=exp_per_gene, pathway.id = "04115", species = "hsa", out.suffix = "p53")
pathview(gene.data=exp_per_gene, pathway.id = "04114", species = "hsa", out.suffix = "oocyte_meiosis")

=cut

use strict;
use warnings;

# REQUIREMENTS
my $group1 = 'State0A'; my $group2 = 'meiosis';
my $degs   = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.degs_$group1"."_vs_$group2"."_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt"; # from 10.run_DEG_for_States0to4.R
my $ids    = '/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/Ens105.ens_to_entrez_ids.txt'; # from BioMart: Gene stable ID, NCBI gene (formerly Entrezgene) ID, Gene name
my $fatal  = 0;
if (!(-e($degs)))  { $fatal++; print "ERROR: cannot find $degs\n";  }
if (!(-e($ids)))   { $fatal++; print "ERROR: cannot find $ids\n";   }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.degs_$group1"."_vs_$group2"."_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.forSPIA.txt";
open(OUT,'>',$out_file) or die $!;
print OUT "EntrezID\tlog2FC\n";

# STORE NAME-TO-ENTREZ GENE ID CONVERSION
my %id_lookup = ();
open(IN,$ids) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/\r$//;
	  my @line = split(/\t/,$line);
	  next if ( (!(defined($line[1]))) or (!(defined($line[2]))) );
	  my $ensembl_gene_id = $line[0]; my $entrez_gene_id = $line[1]; my $gene_name = $line[2];
	  $id_lookup{$gene_name} = $entrez_gene_id;
	}
close(IN) or die $!;

# CREATE A FILE FOR INPUT INTO SPIA, CONTAINING ONLY ENTREZ IDs and LOG2 FOLD CHANGES
my %seen = ();
open(IN,$degs) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line); $line =~ s/[\r\n]//;
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $fc = $line[2]; my $pct1 = $line[3]; my $pct2 = $line[4]; my $adj_p = $line[5];
	  my $diff = abs($pct1-$pct2);
	  next if ($adj_p > 0.05); # CHECKPOINT: we require that DEGs have a significant adj. p-value (< 0.05)
	  next if ($pct2 > $pct1); #CHECKPOINT: we require that each gene be more highly expressed in $group1 than $group2
#	  next if ($diff < 0.2); # CHECKPOINT: we require that genes be differentially detected by >20% of cells
	  next if (exists($seen{$gene_name})); # CHECKPOINT: we have printed this gene to the output file already (this is a necessary checkpoint as some gene names have multiple possible Entrez IDs; we are arbitrarily picking one)
	  if (exists($id_lookup{$gene_name}))
		{ next if ($id_lookup{$gene_name} !~ /^\d+$/);
		  print OUT "$id_lookup{$gene_name}\t$fc\n";
		  $seen{$gene_name}++;
		}
	}
close(IN) or die $!;
close(OUT) or die $!;
exit 1;