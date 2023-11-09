library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(sctransform)
library(sceasy)
library(reticulate)
library(clustree)
theme_set(theme_bw())

### INGEST RAW DATA

# Di Persio 2021

SRS6959440<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS6959440/R_SCT/SRS6959440.rds')
SRS6959441<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS6959441/R_SCT/SRS6959441.rds')
SRS6959442<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS6959442/R_SCT/SRS6959442.rds')

## Guo 2018

SRS3065428<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3065428/R_SCT/SRS3065428.rds')
SRS3065429<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3065429/R_SCT/SRS3065429.rds')
SRS3065430<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3065430/R_SCT/SRS3065430.rds')

## Sohni 2019

SRS4181127<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS4181127/R_SCT/SRS4181127.rds')
SRS4181130<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS4181130/R_SCT/SRS4181130.rds')

## Chen 2020

OES026004<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/OES026004/R_SCT/OES026004.rds')

## Zhao 2020

LZ003<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ003/R_SCT/LZ003.rds')
LZ007<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ007/R_SCT/LZ007.rds')
LZ016<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ016/R_SCT/LZ016.rds')

## Shami 2020

SRS5883810<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883810/R_SCT/SRS5883810.rds')
SRS5883811<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883811/R_SCT/SRS5883811.rds')
SRS5883812<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883812/R_SCT/SRS5883812.rds')
SRS5883813<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883813/R_SCT/SRS5883813.rds')
SRS5883814<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883814/R_SCT/SRS5883814.rds')
SRS5883824<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883824/R_SCT/SRS5883824.rds')
SRS5883825<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883825/R_SCT/SRS5883825.rds')
SRS5883826<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883826/R_SCT/SRS5883826.rds')

## Guo 2020

SRS5086063<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086063/R_SCT/SRS5086063.rds')
SRS5086064<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086064/R_SCT/SRS5086064.rds')

## Nie 2022

SRS9921715<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921715/R_SCT/SRS9921715.rds')
SRS9921716<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921716/R_SCT/SRS9921716.rds')
SRS9921717<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921717/R_SCT/SRS9921717.rds')
SRS9921718<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921718/R_SCT/SRS9921718.rds')
SRS9921719<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921719/R_SCT/SRS9921719.rds')
SRS9921720<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921720/R_SCT/SRS9921720.rds')
SRS9921721<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921721/R_SCT/SRS9921721.rds')
SRS9921722<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921722/R_SCT/SRS9921722.rds')
SRS9921723<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921723/R_SCT/SRS9921723.rds')

## Hermann 2018

SRS2823407<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823407/R_SCT/SRS2823407.rds')
SRS2823408<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823408/R_SCT/SRS2823408.rds')
SRS2823409<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823409/R_SCT/SRS2823409.rds')

### END OF INGESTING RAW DATA

# integrate data using commands from https://github.com/satijalab/seurat/issues/1720 and https://github.com/satijalab/seurat/issues/1836 (the latter recommended by Seurat developers)

testis.list=c(SRS6959440,SRS6959441,SRS6959442,SRS3065428,SRS3065429,SRS3065430,SRS4181127,SRS4181130,OES026004,LZ003,LZ007,LZ016,SRS5883810,SRS5883811,SRS5883812,SRS5883813,SRS5883814,SRS5883824,SRS5883825,SRS5883826,SRS5086063,SRS5086064,SRS9921715,SRS9921716,SRS9921717,SRS9921718,SRS9921719,SRS9921720,SRS9921721,SRS9921722,SRS9921723,SRS2823407,SRS2823408,SRS2823409)
features <- SelectIntegrationFeatures(object.list = testis.list, nfeatures = 5000)
testis.list <- PrepSCTIntegration(object.list = testis.list, anchor.features = features)
testis.anchors <- FindIntegrationAnchors(object.list = testis.list, anchor.features = features, normalization.method = "SCT", reference=c(1,2,3,7,8,10,11,12), reduction = "rpca") # we create as our standard reference the Di Persio, Sohni and Zhao samples as these contain the greatest number of cells
testis.combined <- IntegrateData(anchorset = testis.anchors, normalization.method = "SCT") # see https://github.com/satijalab/seurat/issues/3930 for discussion of k.weight

# add metadata to our integrated data object: cell barcodes, sample IDs, study of origin, ages

testis.combined@meta.data$barcode <- c(
SRS6959440@meta.data$barcode,
SRS6959441@meta.data$barcode,
SRS6959442@meta.data$barcode,
SRS3065428@meta.data$barcode,
SRS3065429@meta.data$barcode,
SRS3065430@meta.data$barcode,
SRS4181127@meta.data$barcode,
SRS4181130@meta.data$barcode,
OES026004@meta.data$barcode,
LZ003@meta.data$barcode,
LZ007@meta.data$barcode,
LZ016@meta.data$barcode,
SRS5883810@meta.data$barcode,
SRS5883811@meta.data$barcode,
SRS5883812@meta.data$barcode,
SRS5883813@meta.data$barcode,
SRS5883814@meta.data$barcode,
SRS5883824@meta.data$barcode,
SRS5883825@meta.data$barcode,
SRS5883826@meta.data$barcode,
SRS5086063@meta.data$barcode,
SRS5086064@meta.data$barcode,
SRS9921715@meta.data$barcode,
SRS9921716@meta.data$barcode,
SRS9921717@meta.data$barcode,
SRS9921718@meta.data$barcode,
SRS9921719@meta.data$barcode,
SRS9921720@meta.data$barcode,
SRS9921721@meta.data$barcode,
SRS9921722@meta.data$barcode,
SRS9921723@meta.data$barcode,
SRS2823407@meta.data$barcode,
SRS2823408@meta.data$barcode,
SRS2823409@meta.data$barcode
)
testis.combined@meta.data$percent.mt <- c(
SRS6959440@meta.data$percent.mt,
SRS6959441@meta.data$percent.mt,
SRS6959442@meta.data$percent.mt,
SRS3065428@meta.data$percent.mt,
SRS3065429@meta.data$percent.mt,
SRS3065430@meta.data$percent.mt,
SRS4181127@meta.data$percent.mt,
SRS4181130@meta.data$percent.mt,
OES026004@meta.data$percent.mt,
LZ003@meta.data$percent.mt,
LZ007@meta.data$percent.mt,
LZ016@meta.data$percent.mt,
SRS5883810@meta.data$percent.mt,
SRS5883811@meta.data$percent.mt,
SRS5883812@meta.data$percent.mt,
SRS5883813@meta.data$percent.mt,
SRS5883814@meta.data$percent.mt,
SRS5883824@meta.data$percent.mt,
SRS5883825@meta.data$percent.mt,
SRS5883826@meta.data$percent.mt,
SRS5086063@meta.data$percent.mt,
SRS5086064@meta.data$percent.mt,
SRS9921715@meta.data$percent.mt,
SRS9921716@meta.data$percent.mt,
SRS9921717@meta.data$percent.mt,
SRS9921718@meta.data$percent.mt,
SRS9921719@meta.data$percent.mt,
SRS9921720@meta.data$percent.mt,
SRS9921721@meta.data$percent.mt,
SRS9921722@meta.data$percent.mt,
SRS9921723@meta.data$percent.mt,
SRS2823407@meta.data$percent.mt,
SRS2823408@meta.data$percent.mt,
SRS2823409@meta.data$percent.mt
)
testis.combined@meta.data$sample.id <- c(
rep("SRS6959440", ncol(SRS6959440@assays$RNA@counts)),
rep("SRS6959441", ncol(SRS6959441@assays$RNA@counts)),
rep("SRS6959442", ncol(SRS6959442@assays$RNA@counts)),
rep("SRS3065428", ncol(SRS3065428@assays$RNA@counts)),
rep("SRS3065429", ncol(SRS3065429@assays$RNA@counts)),
rep("SRS3065430", ncol(SRS3065430@assays$RNA@counts)),
rep("SRS4181127", ncol(SRS4181127@assays$RNA@counts)),
rep("SRS4181130", ncol(SRS4181130@assays$RNA@counts)),
rep("OES026004", ncol(OES026004@assays$RNA@counts)),
rep("LZ003", ncol(LZ003@assays$RNA@counts)),
rep("LZ007", ncol(LZ007@assays$RNA@counts)),
rep("LZ016", ncol(LZ016@assays$RNA@counts)),
rep("SRS5883810", ncol(SRS5883810@assays$RNA@counts)),
rep("SRS5883811", ncol(SRS5883811@assays$RNA@counts)),
rep("SRS5883812", ncol(SRS5883812@assays$RNA@counts)),
rep("SRS5883813", ncol(SRS5883813@assays$RNA@counts)),
rep("SRS5883814", ncol(SRS5883814@assays$RNA@counts)),
rep("SRS5883824", ncol(SRS5883824@assays$RNA@counts)),
rep("SRS5883825", ncol(SRS5883825@assays$RNA@counts)),
rep("SRS5883826", ncol(SRS5883826@assays$RNA@counts)),
rep("SRS5086063", ncol(SRS5086063@assays$RNA@counts)),
rep("SRS5086064", ncol(SRS5086064@assays$RNA@counts)),
rep("SRS9921715", ncol(SRS9921715@assays$RNA@counts)),
rep("SRS9921716", ncol(SRS9921716@assays$RNA@counts)),
rep("SRS9921717", ncol(SRS9921717@assays$RNA@counts)),
rep("SRS9921718", ncol(SRS9921718@assays$RNA@counts)),
rep("SRS9921719", ncol(SRS9921719@assays$RNA@counts)),
rep("SRS9921720", ncol(SRS9921720@assays$RNA@counts)),
rep("SRS9921721", ncol(SRS9921721@assays$RNA@counts)),
rep("SRS9921722", ncol(SRS9921722@assays$RNA@counts)),
rep("SRS9921723", ncol(SRS9921723@assays$RNA@counts)),
rep("SRS2823407", ncol(SRS2823407@assays$RNA@counts)),
rep("SRS2823408", ncol(SRS2823408@assays$RNA@counts)),
rep("SRS2823409", ncol(SRS2823409@assays$RNA@counts))
)
testis.combined@meta.data$source <- c(
rep("Di Persio 2021", ncol(SRS6959440@assays$RNA@counts)),
rep("Di Persio 2021", ncol(SRS6959441@assays$RNA@counts)),
rep("Di Persio 2021", ncol(SRS6959442@assays$RNA@counts)),
rep("Guo 2018", ncol(SRS3065428@assays$RNA@counts)),
rep("Guo 2018", ncol(SRS3065429@assays$RNA@counts)),
rep("Guo 2018", ncol(SRS3065430@assays$RNA@counts)),
rep("Sohni 2019", ncol(SRS4181127@assays$RNA@counts)),
rep("Sohni 2019", ncol(SRS4181130@assays$RNA@counts)),
rep("Chen 2020", ncol(OES026004@assays$RNA@counts)),
rep("Zhao 2020", ncol(LZ003@assays$RNA@counts)),
rep("Zhao 2020", ncol(LZ007@assays$RNA@counts)),
rep("Zhao 2020", ncol(LZ016@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883810@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883811@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883812@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883813@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883814@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883824@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883825@assays$RNA@counts)),
rep("Shami 2020", ncol(SRS5883826@assays$RNA@counts)),
rep("Guo 2020", ncol(SRS5086063@assays$RNA@counts)),
rep("Guo 2020", ncol(SRS5086064@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921715@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921716@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921717@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921718@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921719@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921720@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921721@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921722@assays$RNA@counts)),
rep("Nie 2022", ncol(SRS9921723@assays$RNA@counts)),
rep("Hermann 2018", ncol(SRS2823407@assays$RNA@counts)),
rep("Hermann 2018", ncol(SRS2823408@assays$RNA@counts)),
rep("Hermann 2018", ncol(SRS2823409@assays$RNA@counts))
)
testis.combined@meta.data$age <- c(
rep("31 years", ncol(SRS6959440@assays$RNA@counts)),
rep("33 years", ncol(SRS6959441@assays$RNA@counts)),
rep("55 years", ncol(SRS6959442@assays$RNA@counts)),
rep("17 years", ncol(SRS3065428@assays$RNA@counts)),
rep("24 years", ncol(SRS3065429@assays$RNA@counts)),
rep("25 years", ncol(SRS3065430@assays$RNA@counts)),
rep("37 years", ncol(SRS4181127@assays$RNA@counts)),
rep("42 years", ncol(SRS4181130@assays$RNA@counts)),
rep("36 years", ncol(OES026004@assays$RNA@counts)),
rep("31 years", ncol(LZ003@assays$RNA@counts)),
rep("31 years", ncol(LZ007@assays$RNA@counts)),
rep("17 years", ncol(LZ016@assays$RNA@counts)),
rep("37 years", ncol(SRS5883810@assays$RNA@counts)),
rep("37 years", ncol(SRS5883811@assays$RNA@counts)),
rep("30-40 years", ncol(SRS5883812@assays$RNA@counts)),
rep("30-40 years", ncol(SRS5883813@assays$RNA@counts)),
rep("30-40 years", ncol(SRS5883814@assays$RNA@counts)),
rep("20-25 years", ncol(SRS5883824@assays$RNA@counts)),
rep("20-25 years", ncol(SRS5883825@assays$RNA@counts)),
rep("20-25 years", ncol(SRS5883826@assays$RNA@counts)),
rep("14 years", ncol(SRS5086063@assays$RNA@counts)),
rep("14 years", ncol(SRS5086064@assays$RNA@counts)),
rep("17 years", ncol(SRS9921715@assays$RNA@counts)),
rep("22 years", ncol(SRS9921716@assays$RNA@counts)),
rep("22 years", ncol(SRS9921717@assays$RNA@counts)),
rep("21 years", ncol(SRS9921718@assays$RNA@counts)),
rep("66 years", ncol(SRS9921719@assays$RNA@counts)),
rep("66 years", ncol(SRS9921720@assays$RNA@counts)),
rep("66 years", ncol(SRS9921721@assays$RNA@counts)),
rep("62 years", ncol(SRS9921722@assays$RNA@counts)),
rep("66 years", ncol(SRS9921723@assays$RNA@counts)),
rep("34 years", ncol(SRS2823407@assays$RNA@counts)),
rep("39 years", ncol(SRS2823408@assays$RNA@counts)),
rep("36 years", ncol(SRS2823409@assays$RNA@counts))
)
testis.combined@meta.data$age <- factor(testis.combined@meta.data$age, levels = c("14 years","17 years","20-25 years","21 years","22 years","24 years","25 years","30-40 years","31 years","33 years","34 years","36 years","37 years","39 years","42 years","55 years","62 years","66 years"))

# create a UMAP plot for the combined dataset, part 1: determine dimensionality from the elbow plot

DefaultAssay(testis.combined) <- "integrated"
testis.combined <- RunPCA(testis.combined, verbose = FALSE)

pdf('human_adult_firstPass.elbow_plot.pdf')
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(testis.combined, ndims = 50)
dev.off()

# quantify content of the elbow plot. implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

pct <- testis.combined[["pca"]]@stdev / sum(testis.combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let's take the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data

prin_comp <- min(component1, component2)
write.table(prin_comp,file='human_adult_firstPass.elbow_PC.txt',row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

# create a UMAP plot for the combined dataset, part 2: the plot itself
# see https://github.com/satijalab/seurat/issues/3953: "we recommend the default k=20 for most datasets. As a rule of thumb you do not want to have a higher k than the number of cells in your least populated cell type"
# so we'll fix k but vary the resolution range to experiment with clustering. Be mindful of the comments on clustering made by https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4: "without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types."

resolution.range <- seq(from = 0, to = 0.5, by = 0.05)

testis.combined <- FindNeighbors(testis.combined, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
testis.combined <- FindClusters(testis.combined, algorithm=3, resolution = resolution.range, verbose = FALSE)
testis.combined <- RunUMAP(testis.combined, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)

# print the UMAP plots, coloured according to clusters and to potential batch effects
# we expect to see that there are NO batch effects: that the clusters are not being established on the basis of sample ID, source or age

pdf('human_adult_firstPass.umap.pdf')
tplot = DimPlot(testis.combined, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('human_adult_firstPass.grouped_by_sampleID.pdf')
tplot = DimPlot(testis.combined, reduction = "umap", group.by="sample.id")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('human_adult_firstPass.grouped_by_source.pdf')
tplot = DimPlot(testis.combined, reduction = "umap", group.by="source")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

pdf('human_adult_firstPass.grouped_by_age.pdf')
tplot = DimPlot(testis.combined, reduction = "umap", group.by="age")
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

# print the Clustree plot, so that we may later check which cluster resolution is appropriate
# in the next script, we will re-calculate the clusters according to our manually-chosen resolution

pdf('human_adult_firstPass.clustree_plot.pdf',paper="a4r")
clustree_fig<-clustree(testis.combined,prefix="integrated_snn_res.")
print(clustree_fig)
dev.off()

# remove extraneous metadata from the data object before export
# note that we are not, at this time, removing the variable testis.combined@meta.data$seurat_clusters but it is worth commenting on: this variable contains the results from FindClusters but because we've run FindClusters with a range of values, it will be iteratively overwritten and thereby contain data from the last value in resolution_range, i.e. testis.combined@meta.data$seurat_clusters == testis.combined@meta.data$integrated_snn_res.1.2
# in our next script we will re-run the clustering having decided on an optimal resolution value. It is in that script that we will tidy up the data object further, and export our finalised metadata

testis.combined@meta.data$SCT_snn_res.0<-NULL
testis.combined@meta.data$SCT_snn_res.0.2<-NULL
testis.combined@meta.data$SCT_snn_res.0.4<-NULL
testis.combined@meta.data$SCT_snn_res.0.6<-NULL
testis.combined@meta.data$SCT_snn_res.0.8<-NULL
testis.combined@meta.data$SCT_snn_res.1<-NULL
testis.combined@meta.data$SCT_snn_res.1.2<-NULL
testis.combined@meta.data$SCT_snn_res.1.4<-NULL
testis.combined@meta.data$SCT_snn_res.1.6<-NULL
testis.combined@meta.data$SCT_snn_res.1.8<-NULL
testis.combined@meta.data$SCT_snn_res.2<-NULL
testis.combined@meta.data$orig.ident<-NULL
testis.combined@meta.data$old.ident<-NULL

# save the Seurat object for later use
# note that we are doing this now, BEFORE creating the h5ad, because we are going to scale the dataset (specifically, the 'RNA' data) to do that
# we DON'T want to do that for our SCT/integrated data: "when running sctransform-based workflows, including integration, do not run the ScaleData() function" (see https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1)

saveRDS(testis.combined, file = "human_adult_firstPass.rds")

# convert the Seurat object to an h5ad object for visualisation with cellxgene

DefaultAssay(testis.combined) <- 'RNA'
testis.combined <- NormalizeData(testis.combined)
testis.combined <- FindVariableFeatures(testis.combined, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(testis.combined)
testis <- ScaleData(testis.combined, features = all.genes)
sceasy::convertFormat(testis.combined, from="seurat", to="anndata", outFile='human_adult_firstPass.h5ad')