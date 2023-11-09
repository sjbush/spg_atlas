# BEFORE RUNNING THIS SCRIPT:
# we will have manually reviewed the Clustree plot (human_adult_firstPass.clustree_plot.pdf) and, using cellxgene with the h5ad (human_adult_firstPass.h5ad), the distribution of clusters at each resolution
# commands for doing this: cellxgene launch human_adult_firstPass.h5ad --embedding umap --title whole-testis --disable-diffexp --gene-sets-file biomarkers.csv
# on the basis of this information, we've chosen a resolution optimal for the dataset (0.15), one that neither over- nor underclusters and allows the extraction of the full SSC --> SPG --> primary spermatocyte trajectory
# WHAT THIS SCRIPT DOES:
# this script ingests the 'first pass' integrated dataset, to which a range of cluster resolutions had been applied
# we will now re-run the clustering using our single, empirically derived, resolution value
# we export both a final UMAP and metadata file (supplanting those of the 'first pass' analysis) as well as a table of potential cluster markers (the output of the FindAllMarkers function)
# WHAT THE NEXT STEPS ARE:
# on the basis of prior knowledge and enrichment analysis of the cluster markers, we'll annotate each cluster

library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(sctransform)
library(sceasy)
library(reticulate)
theme_set(theme_bw())

testis.combined<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/human_adult_firstPass.rds')

# purge dataset of existing FindNeighbors/FindClusters data, which is now obsolete

testis.combined@meta.data$integrated_snn_res.0<-NULL
testis.combined@meta.data$integrated_snn_res.0.05<-NULL
testis.combined@meta.data$integrated_snn_res.0.1<-NULL
testis.combined@meta.data$integrated_snn_res.0.15<-NULL
testis.combined@meta.data$integrated_snn_res.0.2<-NULL
testis.combined@meta.data$integrated_snn_res.0.25<-NULL
testis.combined@meta.data$integrated_snn_res.0.3<-NULL
testis.combined@meta.data$integrated_snn_res.0.35<-NULL
testis.combined@meta.data$integrated_snn_res.0.4<-NULL
testis.combined@meta.data$integrated_snn_res.0.45<-NULL
testis.combined@meta.data$integrated_snn_res.0.5<-NULL
testis.combined@meta.data$seurat_clusters<-NULL

# confirm that we are working with the integrated dataset

DefaultAssay(testis.combined) <- "integrated"

# identify the point at which the PCs cover the majority of the variation in the data. We will use this value in the FindNeighors and RunUMAP functions

pct <- testis.combined[["pca"]]@stdev / sum(testis.combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%
prin_comp <- min(component1, component2)

# derive UMAP and clusters

testis.combined <- FindNeighbors(testis.combined, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
testis.combined <- FindClusters(testis.combined, algorithm=3, resolution = 0.15, verbose = FALSE)
testis.combined <- RunUMAP(testis.combined, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)

pdf('/project/GorielyLab2021/sbush/ssc_atlas/human_adult.umap.pdf')
tplot = DimPlot(testis.combined, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# print metadata
# note that testis.combined@meta.data$seurat_clusters == testis.combined@meta.data$integrated_snn_res.0.15, so we don't need to export the latter

testis.combined@meta.data$integrated_snn_res.0.15<-NULL
write.table(testis.combined@meta.data,file='human_adult.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

# determine cluster markers
# we will later run enrichment analyses on these gene lists to assign a predicted function to each cluster
# IMPORTANT: we must set the DefaultAssay to RNA before running differential expression: "we don't recommend using the integrated matrix for differential expression" and "As a general rule, we always recommend performing DE on originally measured values - instead of on batch-corrected, imputed, etc. values. This ensures that the measurements that enter the DE test are indeed independent from each other, which is a requirement of any statistical DE test." (https://github.com/satijalab/seurat/issues/1057, https://github.com/satijalab/seurat/issues/1256 and https://github.com/satijalab/seurat/issues/2136)

DefaultAssay(testis.combined) <- 'RNA'

testis.markers <- FindAllMarkers(testis.combined, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(testis.markers,file='human_adult.cluster_markers.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

# save the Seurat object for later use

saveRDS(testis.combined, file = "human_adult.rds")

# convert the Seurat object to an h5ad object for visualisation with cellxgene

testis.combined <- NormalizeData(testis.combined)
testis.combined <- FindVariableFeatures(testis.combined, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(testis.combined)
testis <- ScaleData(testis.combined, features = all.genes)
sceasy::convertFormat(testis.combined, from="seurat", to="anndata", outFile='human_adult.h5ad')