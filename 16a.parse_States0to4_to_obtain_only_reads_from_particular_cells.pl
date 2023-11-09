=head

PURPOSE:

To calculate RNA velocity, we first need to obtain all reads associated with the cell barcodes used in the hSSC atlas. We only need to calculate velocity for these cells, because we'll later project those velocities into the pre-existing UMAP embedding.
To do this, we export the list of cell barcodes and their associated sample ID from /project/GorielyLab2021/sbush/ssc_atlas/states0to4.rds, the R object created by 8.subset_whole_testes_atlas_to_obtain_States0to4.R

BEFORE USAGE:

module add R-cbrg/202210

library(Seurat)

ssc<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds')

barcode<-rownames(ssc@meta.data)
sample_id<-ssc@meta.data$sample.id
df<-data.frame(BC=barcode,SAMPLE_ID=sample_id)
write.table(df,file='/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

umap <- Embeddings(ssc, reduction = "umap")
write.csv(umap,file='/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.cell_embeddings.csv')

barcode<-rownames(ssc@meta.data)
Idents(ssc)<-ssc$integrated_snn_res.1.1
ssc <- RenameIdents(object = ssc, "0" = "state 0A/1", "1" = "SPG", "2" = "state 0", "3" = "leptotene", "4" = "state 0B", "5" = "zygotene", "6" = "diff SPG")
cluster<-Idents(ssc) # ssc@meta.data$seurat_clusters
df<-data.frame(BC=barcode,CLUSTER=cluster)
write.csv(df, file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.clusters.res1.1.csv",row.names=FALSE)

barcode<-rownames(ssc@meta.data)
Idents(ssc)<-ssc$integrated_snn_res.1.2
ssc <- RenameIdents(object = ssc, "0" = "state 0A/1", "1" = "SPG", "2" = "state 0", "3" = "leptotene", "4" = "zygotene", "5" = "state 0/0B", "6" = "state 0B/1", "7" = "diff SPG")
cluster<-Idents(ssc) # ssc@meta.data$seurat_clusters
df<-data.frame(BC=barcode,CLUSTER=cluster)
write.csv(df, file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.clusters.res1.2.csv",row.names=FALSE)

=cut

use strict;
use warnings;

# REQUIREMENTS
my $in_file = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample.txt"; # created using commands above, from files made by 8.subset_whole_testes_atlas_to_obtain_States0to4.R
if (!(-e($in_file))) { print "ERROR: unable to find $in_file\n"; exit 1; }

# OUTPUT
my $out_dir = "/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample";
if (!(-d($out_dir))) { mkdir $out_dir or die $!; }

my %barcodes_per_sample = ();
open(IN,$in_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $barcode = $line[0]; my $sample_id = $line[1];
	  if ($barcode =~ /^(.+?)\_\d+$/) { $barcode = $1; }
	  $barcodes_per_sample{$sample_id}{$barcode}++;
	}
close(IN) or die $!;
while((my $sample_id,my $irrel)=each(%barcodes_per_sample))
	{ open(OUT,'>',"$out_dir/$sample_id.csv") or die $!;
	  print OUT "\"x\"\n";
	  my @barcodes = ();
	  while((my $barcode,my $irrel)=each(%{$barcodes_per_sample{$sample_id}}))
		  { push(@barcodes,$barcode); }
	  my @sorted_barcodes = sort {$a cmp $b} @barcodes;
	  foreach my $barcode (@sorted_barcodes)
		{ print OUT "\"$barcode\"\n";
		}
	  close(OUT) or die $!;
	}

exit 1;