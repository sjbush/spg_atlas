# AFTER USAGE: visualise using cellxgene launch states0to4.h5ad --embedding umap --title states0to4 --disable-diffexp --max-category-items 500 --gene-sets-file biomarkers.csv

library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(sctransform)
library(reticulate)
library(clustree)

testis<-readRDS('results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
ssc <- subset(testis, (seurat_clusters == 6 | seurat_clusters == 8)) # i.e. retain only cells of states 0-4 (predominantly spermatogonial stem cells to differentiating spermatogonia, although the trajectory ends with some primary spermatocytes)

# IMPORTANT: dataset integration in human_adult.rds is based on variable features in this ORIGINAL, COMPLETE, DATASET
# if we subset that dataset (as we've just done above), we can no longer use that original 'integrated' assay. We will instead have to re-integrate.
# to integrate meaningfully - which requires that we find variable genes - we'll also have to exclude samples with too few cells (< 100). Ultimately, this restricts our dataset to 9 samples
# note that we'll also first have to RE-CALCULATE SCTRANSFORM. Why? From https://github.com/satijalab/seurat/issues/4059: "in the subpopulations, FindVariables will not work for sct-normalized data, because variable genes from sct method is determined by residual variance within the entire population. If you want to re-define variable genes, you may have to re-run SCTransform on your subset populations and integrations."

DefaultAssay(ssc) <- "SCT"
ssc.list <- SplitObject(ssc, split.by = "sample.id")

# exclude 24 samples which contribute < 100 cells to the SSC+SPG supercluster
ssc.list$SRS3065428 <- NULL
ssc.list$SRS3065429 <- NULL
ssc.list$SRS3065430 <- NULL
ssc.list$OES026004 <- NULL
ssc.list$LZ003 <- NULL
ssc.list$LZ016 <- NULL
ssc.list$SRS5883810 <- NULL
ssc.list$SRS5883811 <- NULL
ssc.list$SRS5883812 <- NULL
ssc.list$SRS5883824 <- NULL
ssc.list$SRS5883826 <- NULL
ssc.list$SRS5086063 <- NULL
ssc.list$SRS5086064 <- NULL
ssc.list$SRS9921715 <- NULL
ssc.list$SRS9921716 <- NULL
ssc.list$SRS9921717 <- NULL
ssc.list$SRS9921718 <- NULL
ssc.list$SRS9921719 <- NULL
ssc.list$SRS9921720 <- NULL
ssc.list$SRS9921721 <- NULL
ssc.list$SRS9921722 <- NULL
ssc.list$SRS9921723 <- NULL
ssc.list$SRS2823407 <- NULL

ssc.list <- lapply(X = ssc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ssc.list, nfeatures = 5000)
ssc.list <- PrepSCTIntegration(object.list = ssc.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ssc.list, anchor.features = features, normalization.method = "SCT", reduction = "cca")
full_gene_list<-Reduce(intersect,lapply(ssc.list,rownames))
testis.ssc <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = full_gene_list)

# create a UMAP plot for the State 0-4 dataset, part 1: determine dimensionality from the elbow plot

DefaultAssay(testis.ssc) <- "integrated"
testis.ssc <- RunPCA(testis.ssc, verbose = FALSE)

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

pct <- testis.ssc[["pca"]]@stdev / sum(testis.ssc[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's choose the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data

prin_comp <- min(component1, component2)

# let's initially cluster this data at v. high resolution. We do this because we wish to assign cluster IDs to any (small) contaminating cell populations; this is so we can later remove them

testis.ssc <- FindNeighbors(testis.ssc, reduction = 'pca', dims = 1:prin_comp, k.param = 200, verbose = FALSE)
testis.ssc <- FindClusters(testis.ssc, algorithm=3, resolution = 5.0, verbose = FALSE)
testis.ssc <- RunUMAP(testis.ssc, dims = 1:prin_comp, n.neighbors = 200, verbose = FALSE)

pdf('states0to4_firstPass.umap.resolution_5.0.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# let's save this original, 'first pass', version of the states 0-4 atlas, alongside metadata and clusters
# this is so that we can later scrutinise the contents of any contaminating clusters

saveRDS(testis.ssc, file = "states0to4_firstPass.rds")

write.table(testis.ssc@meta.data,file='states0to4_firstPass.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

#testis.markers <- FindAllMarkers(testis.ssc, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(testis.markers,file='states0to4_firstPass.degs_All_vs_All_resolution_5.0.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

# purge the dataset of the existing FindNeighbors/FindClusters data, as this is now obsolete
# we do this so that when we run saveRDS for the subset, we are not going to export clustering data from the firstPass dataset too
# also: if we didn't remove $integrated_snn_res.5.0, it'd erroneously appear on the Clustree plot too

testis.ssc@meta.data$integrated_snn_res.5<-NULL

# what the above UMAP ("states0to4_firstPass.umap.resolution_5.0.pdf") tells us is that there are some small cell populations (clusters 6, 16 and 17) which are not on the main trajectory - plausible contaminants. We'll remove these and then re-integrate.
# in doing so, we'll need to exclude 2 more samples which no longer contribute > 100 cells to the SSC+SPG supercluster

ssc <- subset(testis.ssc, (seurat_clusters == 6 | seurat_clusters == 16 | seurat_clusters == 17), invert=TRUE)
DefaultAssay(ssc) <- "SCT"
ssc.list <- SplitObject(ssc, split.by = "sample.id")
ssc.list$SRS5883813 <- NULL
ssc.list$SRS5883814 <- NULL
ssc.list <- lapply(X = ssc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ssc.list, nfeatures = 5000)
ssc.list <- PrepSCTIntegration(object.list = ssc.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ssc.list, anchor.features = features, normalization.method = "SCT", reduction = "cca")
full_gene_list<-Reduce(intersect,lapply(ssc.list,rownames))
testis.ssc <- IntegrateData(anchorset = anchors, normalization.method = "SCT", features.to.integrate = full_gene_list)
DefaultAssay(testis.ssc) <- "integrated"
testis.ssc <- RunPCA(testis.ssc, verbose = FALSE)
pct <- testis.ssc[["pca"]]@stdev / sum(testis.ssc[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%
prin_comp <- min(component1, component2)

# we'll now optimise the UMAP by running Clustree across a range of resolutions

resolution.range <- seq(from = 0, to = 2, by = 0.1)

testis.ssc <- FindNeighbors(testis.ssc, reduction = 'pca', dims = 1:prin_comp, k.param = 200, verbose = FALSE)
testis.ssc <- FindClusters(testis.ssc, algorithm=3, resolution = resolution.range, verbose = FALSE)
testis.ssc <- RunUMAP(testis.ssc, dims = 1:prin_comp, n.neighbors = 200, verbose = FALSE)

pdf('states0to4.clustree_plot.pdf',width=7,height=10)
clustree_fig<-clustree(testis.ssc,prefix="integrated_snn_res.")
print(clustree_fig)
dev.off()

# at resolution 1.2, the undifferentiated SSC 'ring' is partitioned into four clusters

testis.ssc <- FindNeighbors(testis.ssc, reduction = 'pca', dims = 1:prin_comp, k.param = 200, verbose = FALSE)
testis.ssc <- FindClusters(testis.ssc, algorithm=3, resolution = 1.2, verbose = FALSE)
testis.ssc <- RunUMAP(testis.ssc, dims = 1:prin_comp, n.neighbors = 200, verbose = FALSE)

pdf('states0to4.umap.resolution_1.2.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# at resolution 1.1, the undifferentiated SSC 'ring' is partitioned into three clusters

testis.ssc <- FindNeighbors(testis.ssc, reduction = 'pca', dims = 1:prin_comp, k.param = 200, verbose = FALSE)
testis.ssc <- FindClusters(testis.ssc, algorithm=3, resolution = 1.1, verbose = FALSE)
testis.ssc <- RunUMAP(testis.ssc, dims = 1:prin_comp, n.neighbors = 200, verbose = FALSE)

pdf('states0to4.umap.resolution_1.1.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# at resolution 0.3, the undifferentiated SSC 'ring' is partitioned into two clusters

testis.ssc <- FindNeighbors(testis.ssc, reduction = 'pca', dims = 1:prin_comp, k.param = 200, verbose = FALSE)
testis.ssc <- FindClusters(testis.ssc, algorithm=3, resolution = 0.3, verbose = FALSE)
testis.ssc <- RunUMAP(testis.ssc, dims = 1:prin_comp, n.neighbors = 200, verbose = FALSE, return.model = TRUE)

pdf('states0to4.umap.resolution_0.3.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('states0to4.grouped_by_sampleID.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('states0to4.grouped_by_source.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", group.by="source")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('states0to4.grouped_by_age.pdf')
tplot = DimPlot(testis.ssc, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

saveRDS(testis.ssc, file = "states0to4.rds")

write.table(testis.ssc@meta.data,file='states0to4.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')