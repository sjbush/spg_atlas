## PROJECTIONS: of a time course of unsorted/unselected mouse samples onto the adult whole testes UMAP. This is a gating strategy to distinguish somatic from germline cells, before re-projecting the latter onto the states 0-4 UMAP
## IMPORTANT: use "module add R-cbrg/202210" to ensure you use v4.1.1 of Seurat

library(Seurat)
library(ggpubr)
library(tidyverse)
theme_set(theme_bw())

# load reference: adult human whole testes data
adult <- readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
adult@meta.data$ssc_state <- Idents(adult)

# to re-create the human UMAP, we first need to quantify the content of the elbow plot, implementing code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- adult[["pca"]]@stdev / sum(adult[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let us choose the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data. We need this parameter to input into RunUMAP
prin_comp <- min(component1, component2)

# load query HRR131888 (embryo, 6 weeks) and then project onto the SSC state 0-4 UMAP
HRR131888<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131888/R_SCT/HRR131888.rds')
DefaultAssay(HRR131888) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131888, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131888 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131888, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131888$predicted.cluster <- factor(HRR131888$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131888 <- DimPlot(HRR131888, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 6 weeks\n(HRR131888)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131888$predicted.cluster))
predicted_sscs.HRR131888 <- subset(HRR131888@meta.data,((HRR131888@meta.data$predicted.cluster != 3) & (HRR131888@meta.data$predicted.cluster != 7) & (HRR131888@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131888,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131888.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131888<-NULL

# load query HRR131892 (embryo, 8 weeks) and then project onto the SSC state 0-4 UMAP
HRR131892<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131892/R_SCT/HRR131892.rds')
DefaultAssay(HRR131892) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131892, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131892 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131892, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131892$predicted.cluster <- factor(HRR131892$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131892 <- DimPlot(HRR131892, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 8 weeks\n(HRR131892)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131892$predicted.cluster))
predicted_sscs.HRR131892 <- subset(HRR131892@meta.data,((HRR131892@meta.data$predicted.cluster != 3) & (HRR131892@meta.data$predicted.cluster != 7) & (HRR131892@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131892,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131892.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131892<-NULL

# load query HRR131894 (embryo, 9 weeks) and then project onto the SSC state 0-4 UMAP
HRR131894<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131894/R_SCT/HRR131894.rds')
DefaultAssay(HRR131894) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131894, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131894 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131894, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131894$predicted.cluster <- factor(HRR131894$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131894 <- DimPlot(HRR131894, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 9 weeks\n(HRR131894)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131894$predicted.cluster))
predicted_sscs.HRR131894 <- subset(HRR131894@meta.data,((HRR131894@meta.data$predicted.cluster != 3) & (HRR131894@meta.data$predicted.cluster != 7) & (HRR131894@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131894,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131894.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131894<-NULL

# load query HRR131895 (embryo, 10 weeks) and then project onto the SSC state 0-4 UMAP
HRR131895<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131895/R_SCT/HRR131895.rds')
DefaultAssay(HRR131895) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131895, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131895 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131895, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131895$predicted.cluster <- factor(HRR131895$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131895 <- DimPlot(HRR131895, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 10 weeks\n(HRR131895)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131895$predicted.cluster))
predicted_sscs.HRR131895 <- subset(HRR131895@meta.data,((HRR131895@meta.data$predicted.cluster != 3) & (HRR131895@meta.data$predicted.cluster != 7) & (HRR131895@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131895,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131895.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131895<-NULL

# load query HRR131897 (embryo, 15 weeks) and then project onto the SSC state 0-4 UMAP
HRR131897<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131897/R_SCT/HRR131897.rds')
DefaultAssay(HRR131897) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131897, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131897 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131897, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131897$predicted.cluster <- factor(HRR131897$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131897 <- DimPlot(HRR131897, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 15 weeks\n(HRR131897)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131897$predicted.cluster))
predicted_sscs.HRR131897 <- subset(HRR131897@meta.data,((HRR131897@meta.data$predicted.cluster != 3) & (HRR131897@meta.data$predicted.cluster != 7) & (HRR131897@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131897,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131897.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131897<-NULL

# load query HRR131901 (embryo, 19 weeks) and then project onto the SSC state 0-4 UMAP
HRR131901<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131901/R_SCT/HRR131901.rds')
DefaultAssay(HRR131901) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131901, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131901 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131901, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131901$predicted.cluster <- factor(HRR131901$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131901 <- DimPlot(HRR131901, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 19 weeks\n(HRR131901)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131901$predicted.cluster))
predicted_sscs.HRR131901 <- subset(HRR131901@meta.data,((HRR131901@meta.data$predicted.cluster != 3) & (HRR131901@meta.data$predicted.cluster != 7) & (HRR131901@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131901,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131901.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131901<-NULL

# load query HRR131904 (embryo, 23 weeks) and then project onto the SSC state 0-4 UMAP
HRR131904<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131904/R_SCT/HRR131904.rds')
DefaultAssay(HRR131904) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131904, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131904 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131904, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131904$predicted.cluster <- factor(HRR131904$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131904 <- DimPlot(HRR131904, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 23 weeks\n(HRR131904)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131904$predicted.cluster))
predicted_sscs.HRR131904 <- subset(HRR131904@meta.data,((HRR131904@meta.data$predicted.cluster != 3) & (HRR131904@meta.data$predicted.cluster != 7) & (HRR131904@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131904,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131904.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131904<-NULL

# load query HRR131906 (embryo, 23 weeks) and then project onto the SSC state 0-4 UMAP
HRR131906<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/HRR131906/R_SCT/HRR131906.rds')
DefaultAssay(HRR131906) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = HRR131906, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
HRR131906 <- MapQuery(anchorset = adult.anchors, reference = adult, query = HRR131906, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
HRR131906$predicted.cluster <- factor(HRR131906$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.HRR131906 <- DimPlot(HRR131906, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("embryo, 23 weeks\n(HRR131906)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(HRR131906$predicted.cluster))
predicted_sscs.HRR131906 <- subset(HRR131906@meta.data,((HRR131906@meta.data$predicted.cluster != 3) & (HRR131906@meta.data$predicted.cluster != 7) & (HRR131906@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.HRR131906,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.HRR131906.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
HRR131906<-NULL

# load query SRS4181123 (2 days) and then project onto the adult whole testes UMAP
SRS4181123<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS4181123/R_SCT/SRS4181123.rds')
DefaultAssay(SRS4181123) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS4181123,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS4181123, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS4181123 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS4181123, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS4181123$predicted.cluster <- factor(SRS4181123$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS4181123 <- DimPlot(SRS4181123, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("2 days\n(SRS4181123)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS4181123$predicted.cluster))
predicted_sscs.SRS4181123 <- subset(SRS4181123@meta.data,((SRS4181123@meta.data$predicted.cluster != 3) & (SRS4181123@meta.data$predicted.cluster != 7) & (SRS4181123@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS4181123,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS4181123.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS4181123<-NULL

# load query SRS4181126 (7 days) and then project onto the adult whole testes UMAP
SRS4181126<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS4181126/R_SCT/SRS4181126.rds')
DefaultAssay(SRS4181126) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS4181126,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS4181126, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS4181126 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS4181126, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS4181126$predicted.cluster <- factor(SRS4181126$predicted.cluster, levels = levels.set)
fig.SRS4181126 <- DimPlot(SRS4181126, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("7 days\n(SRS4181126)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS4181126$predicted.cluster))
predicted_sscs.SRS4181126 <- subset(SRS4181126@meta.data,((SRS4181126@meta.data$predicted.cluster != 3) & (SRS4181126@meta.data$predicted.cluster != 7) & (SRS4181126@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS4181126,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS4181126.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS4181126<-NULL

# load query SRS7727466 (5 months) and then project onto the adult whole testes UMAP
SRS7727466<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS7727466/R_SCT/SRS7727466.rds')
DefaultAssay(SRS7727466) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS7727466,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS7727466, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS7727466 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS7727466, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS7727466$predicted.cluster <- factor(SRS7727466$predicted.cluster, levels = levels.set)
fig.SRS7727466 <- DimPlot(SRS7727466, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 months\n(SRS7727466)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS7727466$predicted.cluster))
predicted_sscs.SRS7727466 <- subset(SRS7727466@meta.data,((SRS7727466@meta.data$predicted.cluster != 3) & (SRS7727466@meta.data$predicted.cluster != 7) & (SRS7727466@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS7727466,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS7727466.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS7727466<-NULL

# load query SRS7727467 (5 months) and then project onto the adult whole testes UMAP
SRS7727467<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS7727467/R_SCT/SRS7727467.rds')
DefaultAssay(SRS7727467) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS7727467,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS7727467, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS7727467 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS7727467, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS7727467$predicted.cluster <- factor(SRS7727467$predicted.cluster, levels = levels.set)
fig.SRS7727467 <- DimPlot(SRS7727467, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 months\n(SRS7727467)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS7727467$predicted.cluster))
predicted_sscs.SRS7727467 <- subset(SRS7727467@meta.data,((SRS7727467@meta.data$predicted.cluster != 3) & (SRS7727467@meta.data$predicted.cluster != 7) & (SRS7727467@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS7727467,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS7727467.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS7727467<-NULL

# load query SRS12015709 (1 year) and then project onto the adult whole testes UMAP
SRS12015709<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS12015709/R_SCT/SRS12015709.rds')
DefaultAssay(SRS12015709) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS12015709,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS12015709, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS12015709 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS12015709, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS12015709$predicted.cluster <- factor(SRS12015709$predicted.cluster, levels = levels.set)
fig.SRS12015709 <- DimPlot(SRS12015709, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("1 year\n(SRS12015709)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS12015709$predicted.cluster))
predicted_sscs.SRS12015709 <- subset(SRS12015709@meta.data,((SRS12015709@meta.data$predicted.cluster != 3) & (SRS12015709@meta.data$predicted.cluster != 7) & (SRS12015709@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS12015709,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS12015709.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS12015709<-NULL

# load query SRS3822680 (13 months) and then project onto the adult whole testes UMAP
SRS3822680<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3822680/R_SCT/SRS3822680.rds')
DefaultAssay(SRS3822680) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3822680,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3822680, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3822680 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3822680, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3822680$predicted.cluster <- factor(SRS3822680$predicted.cluster, levels = levels.set)
fig.SRS3822680 <- DimPlot(SRS3822680, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("13 months\n(SRS3822680)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3822680$predicted.cluster))
predicted_sscs.SRS3822680 <- subset(SRS3822680@meta.data,((SRS3822680@meta.data$predicted.cluster != 3) & (SRS3822680@meta.data$predicted.cluster != 7) & (SRS3822680@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3822680,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3822680.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3822680<-NULL

# load query SRS3822682 (13 months) and then project onto the adult whole testes UMAP
SRS3822682<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3822682/R_SCT/SRS3822682.rds')
DefaultAssay(SRS3822682) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3822682,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3822682, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3822682 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3822682, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3822682$predicted.cluster <- factor(SRS3822682$predicted.cluster, levels = levels.set)
fig.SRS3822682 <- DimPlot(SRS3822682, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("13 months\n(SRS3822682)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3822682$predicted.cluster))
predicted_sscs.SRS3822682 <- subset(SRS3822682@meta.data,((SRS3822682@meta.data$predicted.cluster != 3) & (SRS3822682@meta.data$predicted.cluster != 7) & (SRS3822682@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3822682,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3822682.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3822682<-NULL

# load query SRS3822683 (13 months) and then project onto the adult whole testes UMAP
SRS3822683<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3822683/R_SCT/SRS3822683.rds')
DefaultAssay(SRS3822683) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3822683,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3822683, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3822683 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3822683, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3822683$predicted.cluster <- factor(SRS3822683$predicted.cluster, levels = levels.set)
fig.SRS3822683 <- DimPlot(SRS3822683, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("13 months\n(SRS3822683)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3822683$predicted.cluster))
predicted_sscs.SRS3822683 <- subset(SRS3822683@meta.data,((SRS3822683@meta.data$predicted.cluster != 3) & (SRS3822683@meta.data$predicted.cluster != 7) & (SRS3822683@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3822683,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3822683.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3822683<-NULL

# load query SRS3822686 (13 months) and then project onto the adult whole testes UMAP
SRS3822686<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3822686/R_SCT/SRS3822686.rds')
DefaultAssay(SRS3822686) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3822686,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3822686, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3822686 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3822686, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3822686$predicted.cluster <- factor(SRS3822686$predicted.cluster, levels = levels.set)
fig.SRS3822686 <- DimPlot(SRS3822686, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("13 months\n(SRS3822686)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3822686$predicted.cluster))
predicted_sscs.SRS3822686 <- subset(SRS3822686@meta.data,((SRS3822686@meta.data$predicted.cluster != 3) & (SRS3822686@meta.data$predicted.cluster != 7) & (SRS3822686@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3822686,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3822686.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3822686<-NULL

# load query LZ011 (2 years) and then project onto the adult whole testes UMAP
LZ011<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ011/R_SCT/LZ011.rds')
DefaultAssay(LZ011) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ011,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ011, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ011 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ011, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ011$predicted.cluster <- factor(LZ011$predicted.cluster, levels = levels.set)
fig.LZ011 <- DimPlot(LZ011, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("2 years\n(LZ011)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ011$predicted.cluster))
predicted_sscs.LZ011 <- subset(LZ011@meta.data,((LZ011@meta.data$predicted.cluster != 3) & (LZ011@meta.data$predicted.cluster != 7) & (LZ011@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ011,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ011.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ011<-NULL

# load query SRS12015710 (2 years) and then project onto the adult whole testes UMAP
SRS12015710<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS12015710/R_SCT/SRS12015710.rds')
DefaultAssay(SRS12015710) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS12015710,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS12015710, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS12015710 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS12015710, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS12015710$predicted.cluster <- factor(SRS12015710$predicted.cluster, levels = levels.set)
fig.SRS12015710 <- DimPlot(SRS12015710, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("2 years\n(SRS12015710)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS12015710$predicted.cluster))
predicted_sscs.SRS12015710 <- subset(SRS12015710@meta.data,((SRS12015710@meta.data$predicted.cluster != 3) & (SRS12015710@meta.data$predicted.cluster != 7) & (SRS12015710@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS12015710,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS12015710.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS12015710<-NULL

# load query LZ009 (5 years) and then project onto the adult whole testes UMAP
LZ009<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ009/R_SCT/LZ009.rds')
DefaultAssay(LZ009) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ009,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ009, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ009 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ009, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ009$predicted.cluster <- factor(LZ009$predicted.cluster, levels = levels.set)
fig.LZ009 <- DimPlot(LZ009, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 years\n(LZ009)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ009$predicted.cluster))
predicted_sscs.LZ009 <- subset(LZ009@meta.data,((LZ009@meta.data$predicted.cluster != 3) & (LZ009@meta.data$predicted.cluster != 7) & (LZ009@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ009,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ009.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ009<-NULL

# load query SRS12015708 (7 years) and then project onto the adult whole testes UMAP
SRS12015708<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS12015708/R_SCT/SRS12015708.rds')
DefaultAssay(SRS12015708) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS12015708,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS12015708, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS12015708 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS12015708, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS12015708$predicted.cluster <- factor(SRS12015708$predicted.cluster, levels = levels.set)
fig.SRS12015708 <- DimPlot(SRS12015708, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("7 years\n(SRS12015708)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS12015708$predicted.cluster))
predicted_sscs.SRS12015708 <- subset(SRS12015708@meta.data,((SRS12015708@meta.data$predicted.cluster != 3) & (SRS12015708@meta.data$predicted.cluster != 7) & (SRS12015708@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS12015708,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS12015708.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS12015708<-NULL

# load query SRS5086057 (7 years) and then project onto the adult whole testes UMAP
SRS5086057<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086057/R_SCT/SRS5086057.rds')
DefaultAssay(SRS5086057) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086057,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086057, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086057 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086057, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086057$predicted.cluster <- factor(SRS5086057$predicted.cluster, levels = levels.set)
fig.SRS5086057 <- DimPlot(SRS5086057, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("7 years\n(SRS5086057)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086057$predicted.cluster))
predicted_sscs.SRS5086057 <- subset(SRS5086057@meta.data,((SRS5086057@meta.data$predicted.cluster != 3) & (SRS5086057@meta.data$predicted.cluster != 7) & (SRS5086057@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086057,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086057.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086057<-NULL

# load query SRS5086058 (7 years) and then project onto the adult whole testes UMAP
SRS5086058<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086058/R_SCT/SRS5086058.rds')
DefaultAssay(SRS5086058) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086058,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086058, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086058 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086058, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086058$predicted.cluster <- factor(SRS5086058$predicted.cluster, levels = levels.set)
fig.SRS5086058 <- DimPlot(SRS5086058, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("7 years\n(SRS5086058)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086058$predicted.cluster))
predicted_sscs.SRS5086058 <- subset(SRS5086058@meta.data,((SRS5086058@meta.data$predicted.cluster != 3) & (SRS5086058@meta.data$predicted.cluster != 7) & (SRS5086058@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086058,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086058.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086058<-NULL

# load query LZ005 (8 years) and then project onto the adult whole testes UMAP
LZ005<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ005/R_SCT/LZ005.rds')
DefaultAssay(LZ005) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ005,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ005, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ005 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ005, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ005$predicted.cluster <- factor(LZ005$predicted.cluster, levels = levels.set)
fig.LZ005 <- DimPlot(LZ005, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8 years\n(LZ005)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ005$predicted.cluster))
predicted_sscs.LZ005 <- subset(LZ005@meta.data,((LZ005@meta.data$predicted.cluster != 3) & (LZ005@meta.data$predicted.cluster != 7) & (LZ005@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ005,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ005.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ005<-NULL

# load query LZ008 (11 years) and then project onto the adult whole testes UMAP
LZ008<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ008/R_SCT/LZ008.rds')
DefaultAssay(LZ008) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ008,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ008, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ008 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ008, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ008$predicted.cluster <- factor(LZ008$predicted.cluster, levels = levels.set)
fig.LZ008 <- DimPlot(LZ008, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("11 years\n(LZ008)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ008$predicted.cluster))
predicted_sscs.LZ008 <- subset(LZ008@meta.data,((LZ008@meta.data$predicted.cluster != 3) & (LZ008@meta.data$predicted.cluster != 7) & (LZ008@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ008,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ008.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ008<-NULL

# load query SRS5086059 (11 years) and then project onto the adult whole testes UMAP
SRS5086059<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086059/R_SCT/SRS5086059.rds')
DefaultAssay(SRS5086059) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086059,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086059, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086059 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086059, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086059$predicted.cluster <- factor(SRS5086059$predicted.cluster, levels = levels.set)
fig.SRS5086059 <- DimPlot(SRS5086059, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("11 years\n(SRS5086059)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086059$predicted.cluster))
predicted_sscs.SRS5086059 <- subset(SRS5086059@meta.data,((SRS5086059@meta.data$predicted.cluster != 3) & (SRS5086059@meta.data$predicted.cluster != 7) & (SRS5086059@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086059,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086059.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086059<-NULL

# load query SRS5086060 (11 years) and then project onto the adult whole testes UMAP
SRS5086060<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086060/R_SCT/SRS5086060.rds')
DefaultAssay(SRS5086060) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086060,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086060, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086060 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086060, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086060$predicted.cluster <- factor(SRS5086060$predicted.cluster, levels = levels.set)
fig.SRS5086060 <- DimPlot(SRS5086060, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("11 years\n(SRS5086060)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086060$predicted.cluster))
predicted_sscs.SRS5086060 <- subset(SRS5086060@meta.data,((SRS5086060@meta.data$predicted.cluster != 3) & (SRS5086060@meta.data$predicted.cluster != 7) & (SRS5086060@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086060,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086060.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086060<-NULL

# load query SRS5086061 (13 years) and then project onto the adult whole testes UMAP
SRS5086061<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086061/R_SCT/SRS5086061.rds')
DefaultAssay(SRS5086061) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086061,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086061, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086061 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086061, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086061$predicted.cluster <- factor(SRS5086061$predicted.cluster, levels = levels.set)
fig.SRS5086061 <- DimPlot(SRS5086061, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("13 years\n(SRS5086061)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086061$predicted.cluster))
predicted_sscs.SRS5086061 <- subset(SRS5086061@meta.data,((SRS5086061@meta.data$predicted.cluster != 3) & (SRS5086061@meta.data$predicted.cluster != 7) & (SRS5086061@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086061,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086061.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086061<-NULL

# load query SRS5086062 (13 years) and then project onto the adult whole testes UMAP
SRS5086062<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086062/R_SCT/SRS5086062.rds')
DefaultAssay(SRS5086062) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086062,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086062, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086062 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086062, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086062$predicted.cluster <- factor(SRS5086062$predicted.cluster, levels = levels.set)
fig.SRS5086062 <- DimPlot(SRS5086062, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("13 years\n(SRS5086062)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086062$predicted.cluster))
predicted_sscs.SRS5086062 <- subset(SRS5086062@meta.data,((SRS5086062@meta.data$predicted.cluster != 3) & (SRS5086062@meta.data$predicted.cluster != 7) & (SRS5086062@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086062,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086062.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086062<-NULL

# load query SRS5086063 (14 years) and then project onto the adult whole testes UMAP
SRS5086063<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086063/R_SCT/SRS5086063.rds')
DefaultAssay(SRS5086063) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086063,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086063, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086063 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086063, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086063$predicted.cluster <- factor(SRS5086063$predicted.cluster, levels = levels.set)
fig.SRS5086063 <- DimPlot(SRS5086063, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("14 years\n(SRS5086063)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086063$predicted.cluster))
predicted_sscs.SRS5086063 <- subset(SRS5086063@meta.data,((SRS5086063@meta.data$predicted.cluster != 3) & (SRS5086063@meta.data$predicted.cluster != 7) & (SRS5086063@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086063,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086063.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086063<-NULL

# load query SRS5086064 (14 years) and then project onto the adult whole testes UMAP
SRS5086064<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086064/R_SCT/SRS5086064.rds')
DefaultAssay(SRS5086064) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086064,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086064, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086064 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086064, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086064$predicted.cluster <- factor(SRS5086064$predicted.cluster, levels = levels.set)
fig.SRS5086064 <- DimPlot(SRS5086064, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("14 years\n(SRS5086064)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086064$predicted.cluster))
predicted_sscs.SRS5086064 <- subset(SRS5086064@meta.data,((SRS5086064@meta.data$predicted.cluster != 3) & (SRS5086064@meta.data$predicted.cluster != 7) & (SRS5086064@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086064,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086064.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086064<-NULL

# load query SRS3065428 (17 years) and then project onto the adult whole testes UMAP
SRS3065428<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3065428/R_SCT/SRS3065428.rds')
DefaultAssay(SRS3065428) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3065428,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3065428, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3065428 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3065428, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3065428$predicted.cluster <- factor(SRS3065428$predicted.cluster, levels = levels.set)
fig.SRS3065428 <- DimPlot(SRS3065428, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("17 years\n(SRS3065428)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3065428$predicted.cluster))
predicted_sscs.SRS3065428 <- subset(SRS3065428@meta.data,((SRS3065428@meta.data$predicted.cluster != 3) & (SRS3065428@meta.data$predicted.cluster != 7) & (SRS3065428@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3065428,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3065428.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3065428<-NULL

# load query LZ016 (17 years) and then project onto the adult whole testes UMAP
LZ016<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ016/R_SCT/LZ016.rds')
DefaultAssay(LZ016) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ016,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ016, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ016 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ016, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ016$predicted.cluster <- factor(LZ016$predicted.cluster, levels = levels.set)
fig.LZ016 <- DimPlot(LZ016, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("17 years\n(LZ016)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ016$predicted.cluster))
predicted_sscs.LZ016 <- subset(LZ016@meta.data,((LZ016@meta.data$predicted.cluster != 3) & (LZ016@meta.data$predicted.cluster != 7) & (LZ016@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ016,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ016.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ016<-NULL

# load query SRS9921715 (17 years) and then project onto the adult whole testes UMAP
SRS9921715<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921715/R_SCT/SRS9921715.rds')
DefaultAssay(SRS9921715) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921715,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921715, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921715 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921715, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921715$predicted.cluster <- factor(SRS9921715$predicted.cluster, levels = levels.set)
fig.SRS9921715 <- DimPlot(SRS9921715, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("17 years\n(SRS9921715)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921715$predicted.cluster))
predicted_sscs.SRS9921715 <- subset(SRS9921715@meta.data,((SRS9921715@meta.data$predicted.cluster != 3) & (SRS9921715@meta.data$predicted.cluster != 7) & (SRS9921715@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921715,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921715.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921715<-NULL

# load query SRS5883824 (20-25 years) and then project onto the adult whole testes UMAP
SRS5883824<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883824/R_SCT/SRS5883824.rds')
DefaultAssay(SRS5883824) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883824,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883824, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883824 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883824, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883824$predicted.cluster <- factor(SRS5883824$predicted.cluster, levels = levels.set)
fig.SRS5883824 <- DimPlot(SRS5883824, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("20-25 years\n(SRS5883824)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883824$predicted.cluster))
predicted_sscs.SRS5883824 <- subset(SRS5883824@meta.data,((SRS5883824@meta.data$predicted.cluster != 3) & (SRS5883824@meta.data$predicted.cluster != 7) & (SRS5883824@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883824,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883824.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883824<-NULL

# load query SRS5883825 (20-25 years) and then project onto the adult whole testes UMAP
SRS5883825<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883825/R_SCT/SRS5883825.rds')
DefaultAssay(SRS5883825) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883825,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883825, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883825 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883825, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883825$predicted.cluster <- factor(SRS5883825$predicted.cluster, levels = levels.set)
fig.SRS5883825 <- DimPlot(SRS5883825, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("20-25 years\n(SRS5883825)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883825$predicted.cluster))
predicted_sscs.SRS5883825 <- subset(SRS5883825@meta.data,((SRS5883825@meta.data$predicted.cluster != 3) & (SRS5883825@meta.data$predicted.cluster != 7) & (SRS5883825@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883825,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883825.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883825<-NULL

# load query SRS5883826 (20-25 years) and then project onto the adult whole testes UMAP
SRS5883826<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883826/R_SCT/SRS5883826.rds')
DefaultAssay(SRS5883826) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883826,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883826, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883826 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883826, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883826$predicted.cluster <- factor(SRS5883826$predicted.cluster, levels = levels.set)
fig.SRS5883826 <- DimPlot(SRS5883826, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("20-25 years\n(SRS5883826)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883826$predicted.cluster))
predicted_sscs.SRS5883826 <- subset(SRS5883826@meta.data,((SRS5883826@meta.data$predicted.cluster != 3) & (SRS5883826@meta.data$predicted.cluster != 7) & (SRS5883826@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883826,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883826.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883826<-NULL

# load query SRS9921718 (21 years) and then project onto the adult whole testes UMAP
SRS9921718<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921718/R_SCT/SRS9921718.rds')
DefaultAssay(SRS9921718) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921718,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921718, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921718 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921718, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921718$predicted.cluster <- factor(SRS9921718$predicted.cluster, levels = levels.set)
fig.SRS9921718 <- DimPlot(SRS9921718, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("21 years\n(SRS9921718)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921718$predicted.cluster))
predicted_sscs.SRS9921718 <- subset(SRS9921718@meta.data,((SRS9921718@meta.data$predicted.cluster != 3) & (SRS9921718@meta.data$predicted.cluster != 7) & (SRS9921718@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921718,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921718.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921718<-NULL

# load query SRS9921716 (22 years) and then project onto the adult whole testes UMAP
SRS9921716<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921716/R_SCT/SRS9921716.rds')
DefaultAssay(SRS9921716) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921716,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921716, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921716 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921716, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921716$predicted.cluster <- factor(SRS9921716$predicted.cluster, levels = levels.set)
fig.SRS9921716 <- DimPlot(SRS9921716, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("22 years\n(SRS9921716)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921716$predicted.cluster))
predicted_sscs.SRS9921716 <- subset(SRS9921716@meta.data,((SRS9921716@meta.data$predicted.cluster != 3) & (SRS9921716@meta.data$predicted.cluster != 7) & (SRS9921716@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921716,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921716.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921716<-NULL

# load query SRS9921717 (22 years) and then project onto the adult whole testes UMAP
SRS9921717<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921717/R_SCT/SRS9921717.rds')
DefaultAssay(SRS9921717) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921717,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921717, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921717 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921717, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921717$predicted.cluster <- factor(SRS9921717$predicted.cluster, levels = levels.set)
fig.SRS9921717 <- DimPlot(SRS9921717, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("22 years\n(SRS9921717)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921717$predicted.cluster))
predicted_sscs.SRS9921717 <- subset(SRS9921717@meta.data,((SRS9921717@meta.data$predicted.cluster != 3) & (SRS9921717@meta.data$predicted.cluster != 7) & (SRS9921717@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921717,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921717.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921717<-NULL

# load query SRS3065429 (24 years) and then project onto the adult whole testes UMAP
SRS3065429<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3065429/R_SCT/SRS3065429.rds')
DefaultAssay(SRS3065429) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3065429,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3065429, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3065429 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3065429, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3065429$predicted.cluster <- factor(SRS3065429$predicted.cluster, levels = levels.set)
fig.SRS3065429 <- DimPlot(SRS3065429, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("24 years\n(SRS3065429)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3065429$predicted.cluster))
predicted_sscs.SRS3065429 <- subset(SRS3065429@meta.data,((SRS3065429@meta.data$predicted.cluster != 3) & (SRS3065429@meta.data$predicted.cluster != 7) & (SRS3065429@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3065429,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3065429.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3065429<-NULL

# load query SRS3065430 (25 years) and then project onto the adult whole testes UMAP
SRS3065430<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS3065430/R_SCT/SRS3065430.rds')
DefaultAssay(SRS3065430) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS3065430,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3065430, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3065430 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3065430, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3065430$predicted.cluster <- factor(SRS3065430$predicted.cluster, levels = levels.set)
fig.SRS3065430 <- DimPlot(SRS3065430, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("25 years\n(SRS3065430)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3065430$predicted.cluster))
predicted_sscs.SRS3065430 <- subset(SRS3065430@meta.data,((SRS3065430@meta.data$predicted.cluster != 3) & (SRS3065430@meta.data$predicted.cluster != 7) & (SRS3065430@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3065430,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3065430.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3065430<-NULL

# load query SRS5787712 (26 years - transfemale) and then project onto the adult whole testes UMAP
SRS5787712<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5787712/R_SCT/SRS5787712.rds')
DefaultAssay(SRS5787712) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5787712,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5787712, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5787712 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5787712, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5787712$predicted.cluster <- factor(SRS5787712$predicted.cluster, levels = levels.set)
fig.SRS5787712 <- DimPlot(SRS5787712, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("26 years - transfemale\n(SRS5787712)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5787712$predicted.cluster))
predicted_sscs.SRS5787712 <- subset(SRS5787712@meta.data,((SRS5787712@meta.data$predicted.cluster != 3) & (SRS5787712@meta.data$predicted.cluster != 7) & (SRS5787712@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5787712,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5787712.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5787712<-NULL

# load query SRS5883812 (30-40 years) and then project onto the adult whole testes UMAP
SRS5883812<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883812/R_SCT/SRS5883812.rds')
DefaultAssay(SRS5883812) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883812,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883812, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883812 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883812, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883812$predicted.cluster <- factor(SRS5883812$predicted.cluster, levels = levels.set)
fig.SRS5883812 <- DimPlot(SRS5883812, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("30-40 years\n(SRS5883812)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883812$predicted.cluster))
predicted_sscs.SRS5883812 <- subset(SRS5883812@meta.data,((SRS5883812@meta.data$predicted.cluster != 3) & (SRS5883812@meta.data$predicted.cluster != 7) & (SRS5883812@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883812,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883812.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883812<-NULL

# load query SRS5883813 (30-40 years) and then project onto the adult whole testes UMAP
SRS5883813<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883813/R_SCT/SRS5883813.rds')
DefaultAssay(SRS5883813) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883813,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883813, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883813 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883813, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883813$predicted.cluster <- factor(SRS5883813$predicted.cluster, levels = levels.set)
fig.SRS5883813 <- DimPlot(SRS5883813, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("30-40 years\n(SRS5883813)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883813$predicted.cluster))
predicted_sscs.SRS5883813 <- subset(SRS5883813@meta.data,((SRS5883813@meta.data$predicted.cluster != 3) & (SRS5883813@meta.data$predicted.cluster != 7) & (SRS5883813@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883813,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883813.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883813<-NULL

# load query SRS5883814 (30-40 years) and then project onto the adult whole testes UMAP
SRS5883814<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883814/R_SCT/SRS5883814.rds')
DefaultAssay(SRS5883814) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883814,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883814, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883814 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883814, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883814$predicted.cluster <- factor(SRS5883814$predicted.cluster, levels = levels.set)
fig.SRS5883814 <- DimPlot(SRS5883814, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("30-40 years\n(SRS5883814)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883814$predicted.cluster))
predicted_sscs.SRS5883814 <- subset(SRS5883814@meta.data,((SRS5883814@meta.data$predicted.cluster != 3) & (SRS5883814@meta.data$predicted.cluster != 7) & (SRS5883814@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883814,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883814.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883814<-NULL

# load query LZ003 (31 years) and then project onto the adult whole testes UMAP
LZ003<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ003/R_SCT/LZ003.rds')
DefaultAssay(LZ003) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ003,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ003, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ003 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ003, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ003$predicted.cluster <- factor(LZ003$predicted.cluster, levels = levels.set)
fig.LZ003 <- DimPlot(LZ003, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("31 years\n(LZ003)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ003$predicted.cluster))
predicted_sscs.LZ003 <- subset(LZ003@meta.data,((LZ003@meta.data$predicted.cluster != 3) & (LZ003@meta.data$predicted.cluster != 7) & (LZ003@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ003,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ003.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ003<-NULL

# load query LZ007 (31 years) and then project onto the adult whole testes UMAP
LZ007<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/LZ007/R_SCT/LZ007.rds')
DefaultAssay(LZ007) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(LZ007,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = LZ007, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
LZ007 <- MapQuery(anchorset = adult.anchors, reference = adult, query = LZ007, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
LZ007$predicted.cluster <- factor(LZ007$predicted.cluster, levels = levels.set)
fig.LZ007 <- DimPlot(LZ007, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("31 years\n(LZ007)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(LZ007$predicted.cluster))
predicted_sscs.LZ007 <- subset(LZ007@meta.data,((LZ007@meta.data$predicted.cluster != 3) & (LZ007@meta.data$predicted.cluster != 7) & (LZ007@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.LZ007,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.LZ007.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
LZ007<-NULL

# load query SRS6959440 (31 years) and then project onto the adult whole testes UMAP
SRS6959440<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS6959440/R_SCT/SRS6959440.rds')
DefaultAssay(SRS6959440) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS6959440,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS6959440, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS6959440 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS6959440, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS6959440$predicted.cluster <- factor(SRS6959440$predicted.cluster, levels = levels.set)
fig.SRS6959440 <- DimPlot(SRS6959440, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("31 years\n(SRS6959440)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS6959440$predicted.cluster))
predicted_sscs.SRS6959440 <- subset(SRS6959440@meta.data,((SRS6959440@meta.data$predicted.cluster != 3) & (SRS6959440@meta.data$predicted.cluster != 7) & (SRS6959440@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS6959440,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS6959440.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS6959440<-NULL

# load query SRS6959441 (33 years) and then project onto the adult whole testes UMAP
SRS6959441<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS6959441/R_SCT/SRS6959441.rds')
DefaultAssay(SRS6959441) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS6959441,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS6959441, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS6959441 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS6959441, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS6959441$predicted.cluster <- factor(SRS6959441$predicted.cluster, levels = levels.set)
fig.SRS6959441 <- DimPlot(SRS6959441, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("33 years\n(SRS6959441)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS6959441$predicted.cluster))
predicted_sscs.SRS6959441 <- subset(SRS6959441@meta.data,((SRS6959441@meta.data$predicted.cluster != 3) & (SRS6959441@meta.data$predicted.cluster != 7) & (SRS6959441@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS6959441,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS6959441.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS6959441<-NULL

# load query SRS2823407 (34 years) and then project onto the adult whole testes UMAP
SRS2823407<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823407/R_SCT/SRS2823407.rds')
DefaultAssay(SRS2823407) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS2823407,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS2823407, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS2823407 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS2823407, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS2823407$predicted.cluster <- factor(SRS2823407$predicted.cluster, levels = levels.set)
fig.SRS2823407 <- DimPlot(SRS2823407, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("34 years\n(SRS2823407)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS2823407$predicted.cluster))
predicted_sscs.SRS2823407 <- subset(SRS2823407@meta.data,((SRS2823407@meta.data$predicted.cluster != 3) & (SRS2823407@meta.data$predicted.cluster != 7) & (SRS2823407@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS2823407,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS2823407.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS2823407<-NULL

# load query SRS2823409 (36 years) and then project onto the adult whole testes UMAP
SRS2823409<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823409/R_SCT/SRS2823409.rds')
DefaultAssay(SRS2823409) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS2823409,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS2823409, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS2823409 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS2823409, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS2823409$predicted.cluster <- factor(SRS2823409$predicted.cluster, levels = levels.set)
fig.SRS2823409 <- DimPlot(SRS2823409, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("36 years\n(SRS2823409)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS2823409$predicted.cluster))
predicted_sscs.SRS2823409 <- subset(SRS2823409@meta.data,((SRS2823409@meta.data$predicted.cluster != 3) & (SRS2823409@meta.data$predicted.cluster != 7) & (SRS2823409@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS2823409,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS2823409.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS2823409<-NULL

# load query OES026004 (36 years) and then project onto the adult whole testes UMAP
OES026004<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/OES026004/R_SCT/OES026004.rds')
DefaultAssay(OES026004) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(OES026004,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = OES026004, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
OES026004 <- MapQuery(anchorset = adult.anchors, reference = adult, query = OES026004, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
OES026004$predicted.cluster <- factor(OES026004$predicted.cluster, levels = levels.set)
fig.OES026004 <- DimPlot(OES026004, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("36 years\n(OES026004)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(OES026004$predicted.cluster))
predicted_sscs.OES026004 <- subset(OES026004@meta.data,((OES026004@meta.data$predicted.cluster != 3) & (OES026004@meta.data$predicted.cluster != 7) & (OES026004@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.OES026004,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.OES026004.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
OES026004<-NULL

# load query SRS5883810 (37 years) and then project onto the adult whole testes UMAP
SRS5883810<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883810/R_SCT/SRS5883810.rds')
DefaultAssay(SRS5883810) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883810,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883810, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883810 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883810, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883810$predicted.cluster <- factor(SRS5883810$predicted.cluster, levels = levels.set)
fig.SRS5883810 <- DimPlot(SRS5883810, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("37 years\n(SRS5883810)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883810$predicted.cluster))
predicted_sscs.SRS5883810 <- subset(SRS5883810@meta.data,((SRS5883810@meta.data$predicted.cluster != 3) & (SRS5883810@meta.data$predicted.cluster != 7) & (SRS5883810@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883810,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883810.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883810<-NULL

# load query SRS5883811 (37 years) and then project onto the adult whole testes UMAP
SRS5883811<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5883811/R_SCT/SRS5883811.rds')
DefaultAssay(SRS5883811) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5883811,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5883811, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5883811 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5883811, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5883811$predicted.cluster <- factor(SRS5883811$predicted.cluster, levels = levels.set)
fig.SRS5883811 <- DimPlot(SRS5883811, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("37 years\n(SRS5883811)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5883811$predicted.cluster))
predicted_sscs.SRS5883811 <- subset(SRS5883811@meta.data,((SRS5883811@meta.data$predicted.cluster != 3) & (SRS5883811@meta.data$predicted.cluster != 7) & (SRS5883811@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5883811,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5883811.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5883811<-NULL

# load query SRS4181127 (37 years) and then project onto the adult whole testes UMAP
SRS4181127<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS4181127/R_SCT/SRS4181127.rds')
DefaultAssay(SRS4181127) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS4181127,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS4181127, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS4181127 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS4181127, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS4181127$predicted.cluster <- factor(SRS4181127$predicted.cluster, levels = levels.set)
fig.SRS4181127 <- DimPlot(SRS4181127, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("37 years\n(SRS4181127)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS4181127$predicted.cluster))
predicted_sscs.SRS4181127 <- subset(SRS4181127@meta.data,((SRS4181127@meta.data$predicted.cluster != 3) & (SRS4181127@meta.data$predicted.cluster != 7) & (SRS4181127@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS4181127,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS4181127.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS4181127<-NULL

# load query SRS4181130 (42 years) and then project onto the adult whole testes UMAP
SRS4181130<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS4181130/R_SCT/SRS4181130.rds')
DefaultAssay(SRS4181130) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS4181130,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS4181130, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS4181130 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS4181130, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS4181130$predicted.cluster <- factor(SRS4181130$predicted.cluster, levels = levels.set)
fig.SRS4181130 <- DimPlot(SRS4181130, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("42 years\n(SRS4181130)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS4181130$predicted.cluster))
predicted_sscs.SRS4181130 <- subset(SRS4181130@meta.data,((SRS4181130@meta.data$predicted.cluster != 3) & (SRS4181130@meta.data$predicted.cluster != 7) & (SRS4181130@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS4181130,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS4181130.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS4181130<-NULL

# load query SRS2823408 (49 years) and then project onto the adult whole testes UMAP
SRS2823408<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823408/R_SCT/SRS2823408.rds')
DefaultAssay(SRS2823408) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS2823408,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS2823408, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS2823408 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS2823408, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS2823408$predicted.cluster <- factor(SRS2823408$predicted.cluster, levels = levels.set)
fig.SRS2823408 <- DimPlot(SRS2823408, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("49 years\n(SRS2823408)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS2823408$predicted.cluster))
predicted_sscs.SRS2823408 <- subset(SRS2823408@meta.data,((SRS2823408@meta.data$predicted.cluster != 3) & (SRS2823408@meta.data$predicted.cluster != 7) & (SRS2823408@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS2823408,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS2823408.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS2823408<-NULL

# load query SRS5086065 (50 years - transfemale) and then project onto the adult whole testes UMAP
SRS5086065<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086065/R_SCT/SRS5086065.rds')
DefaultAssay(SRS5086065) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086065,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086065, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086065 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086065, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086065$predicted.cluster <- factor(SRS5086065$predicted.cluster, levels = levels.set)
fig.SRS5086065 <- DimPlot(SRS5086065, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("50 years - transfemale\n(SRS5086065)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086065$predicted.cluster))
predicted_sscs.SRS5086065 <- subset(SRS5086065@meta.data,((SRS5086065@meta.data$predicted.cluster != 3) & (SRS5086065@meta.data$predicted.cluster != 7) & (SRS5086065@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086065,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086065.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086065<-NULL

# load query SRS5086066 (50 years - transfemale) and then project onto the adult whole testes UMAP
SRS5086066<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS5086066/R_SCT/SRS5086066.rds')
DefaultAssay(SRS5086066) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS5086066,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS5086066, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS5086066 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS5086066, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS5086066$predicted.cluster <- factor(SRS5086066$predicted.cluster, levels = levels.set)
fig.SRS5086066 <- DimPlot(SRS5086066, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("50 years - transfemale\n(SRS5086066)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS5086066$predicted.cluster))
predicted_sscs.SRS5086066 <- subset(SRS5086066@meta.data,((SRS5086066@meta.data$predicted.cluster != 3) & (SRS5086066@meta.data$predicted.cluster != 7) & (SRS5086066@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS5086066,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS5086066.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS5086066<-NULL

# load query SRS6959442 (55 years) and then project onto the adult whole testes UMAP
SRS6959442<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS6959442/R_SCT/SRS6959442.rds')
DefaultAssay(SRS6959442) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS6959442,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS6959442, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS6959442 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS6959442, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS6959442$predicted.cluster <- factor(SRS6959442$predicted.cluster, levels = levels.set)
fig.SRS6959442 <- DimPlot(SRS6959442, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("55 years\n(SRS6959442)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS6959442$predicted.cluster))
predicted_sscs.SRS6959442 <- subset(SRS6959442@meta.data,((SRS6959442@meta.data$predicted.cluster != 3) & (SRS6959442@meta.data$predicted.cluster != 7) & (SRS6959442@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS6959442,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS6959442.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS6959442<-NULL

# load query SRS9921722 (62 years) and then project onto the adult whole testes UMAP
SRS9921722<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921722/R_SCT/SRS9921722.rds')
DefaultAssay(SRS9921722) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921722,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921722, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921722 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921722, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921722$predicted.cluster <- factor(SRS9921722$predicted.cluster, levels = levels.set)
fig.SRS9921722 <- DimPlot(SRS9921722, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("62 years\n(SRS9921722)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921722$predicted.cluster))
predicted_sscs.SRS9921722 <- subset(SRS9921722@meta.data,((SRS9921722@meta.data$predicted.cluster != 3) & (SRS9921722@meta.data$predicted.cluster != 7) & (SRS9921722@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921722,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921722.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921722<-NULL

# load query SRS9921724 (64 years - high BMI) and then project onto the adult whole testes UMAP
SRS9921724<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921724/R_SCT/SRS9921724.rds')
DefaultAssay(SRS9921724) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921724,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921724, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921724 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921724, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921724$predicted.cluster <- factor(SRS9921724$predicted.cluster, levels = levels.set)
fig.SRS9921724 <- DimPlot(SRS9921724, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("64 years - high BMI\n(SRS9921724)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921724$predicted.cluster))
predicted_sscs.SRS9921724 <- subset(SRS9921724@meta.data,((SRS9921724@meta.data$predicted.cluster != 3) & (SRS9921724@meta.data$predicted.cluster != 7) & (SRS9921724@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921724,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921724.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921724<-NULL

# load query SRS9921719 (66 years) and then project onto the adult whole testes UMAP
SRS9921719<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921719/R_SCT/SRS9921719.rds')
DefaultAssay(SRS9921719) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921719,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921719, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921719 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921719, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921719$predicted.cluster <- factor(SRS9921719$predicted.cluster, levels = levels.set)
fig.SRS9921719 <- DimPlot(SRS9921719, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("66 years\n(SRS9921719)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921719$predicted.cluster))
predicted_sscs.SRS9921719 <- subset(SRS9921719@meta.data,((SRS9921719@meta.data$predicted.cluster != 3) & (SRS9921719@meta.data$predicted.cluster != 7) & (SRS9921719@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921719,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921719.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921719<-NULL

# load query SRS9921720 (66 years) and then project onto the adult whole testes UMAP
SRS9921720<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921720/R_SCT/SRS9921720.rds')
DefaultAssay(SRS9921720) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921720,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921720, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921720 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921720, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921720$predicted.cluster <- factor(SRS9921720$predicted.cluster, levels = levels.set)
fig.SRS9921720 <- DimPlot(SRS9921720, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("66 years\n(SRS9921720)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921720$predicted.cluster))
predicted_sscs.SRS9921720 <- subset(SRS9921720@meta.data,((SRS9921720@meta.data$predicted.cluster != 3) & (SRS9921720@meta.data$predicted.cluster != 7) & (SRS9921720@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921720,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921720.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921720<-NULL

# load query SRS9921721 (66 years) and then project onto the adult whole testes UMAP
SRS9921721<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921721/R_SCT/SRS9921721.rds')
DefaultAssay(SRS9921721) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921721,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921721, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921721 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921721, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921721$predicted.cluster <- factor(SRS9921721$predicted.cluster, levels = levels.set)
fig.SRS9921721 <- DimPlot(SRS9921721, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("66 years\n(SRS9921721)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921721$predicted.cluster))
predicted_sscs.SRS9921721 <- subset(SRS9921721@meta.data,((SRS9921721@meta.data$predicted.cluster != 3) & (SRS9921721@meta.data$predicted.cluster != 7) & (SRS9921721@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921721,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921721.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921721<-NULL

# load query SRS9921723 (66 years) and then project onto the adult whole testes UMAP
SRS9921723<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921723/R_SCT/SRS9921723.rds')
DefaultAssay(SRS9921723) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921723,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921723, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921723 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921723, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921723$predicted.cluster <- factor(SRS9921723$predicted.cluster, levels = levels.set)
fig.SRS9921723 <- DimPlot(SRS9921723, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("66 years\n(SRS9921723)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921723$predicted.cluster))
predicted_sscs.SRS9921723 <- subset(SRS9921723@meta.data,((SRS9921723@meta.data$predicted.cluster != 3) & (SRS9921723@meta.data$predicted.cluster != 7) & (SRS9921723@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921723,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921723.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921723<-NULL

# load query SRS9921726 (67 years - high BMI) and then project onto the adult whole testes UMAP
SRS9921726<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921726/R_SCT/SRS9921726.rds')
DefaultAssay(SRS9921726) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921726,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921726, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921726 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921726, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921726$predicted.cluster <- factor(SRS9921726$predicted.cluster, levels = levels.set)
fig.SRS9921726 <- DimPlot(SRS9921726, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("67 years - high BMI\n(SRS9921726)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921726$predicted.cluster))
predicted_sscs.SRS9921726 <- subset(SRS9921726@meta.data,((SRS9921726@meta.data$predicted.cluster != 3) & (SRS9921726@meta.data$predicted.cluster != 7) & (SRS9921726@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921726,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921726.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921726<-NULL

# load query SRS9921725 (76 years - high BMI) and then project onto the adult whole testes UMAP
SRS9921725<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS9921725/R_SCT/SRS9921725.rds')
DefaultAssay(SRS9921725) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS9921725,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS9921725, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS9921725 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS9921725, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS9921725$predicted.cluster <- factor(SRS9921725$predicted.cluster, levels = levels.set)
fig.SRS9921725 <- DimPlot(SRS9921725, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("76 years - high BMI\n(SRS9921725)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS9921725$predicted.cluster))
predicted_sscs.SRS9921725 <- subset(SRS9921725@meta.data,((SRS9921725@meta.data$predicted.cluster != 3) & (SRS9921725@meta.data$predicted.cluster != 7) & (SRS9921725@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS9921725,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9921725.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9921725<-NULL

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt1.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.HRR131888,fig.HRR131892,fig.HRR131894,fig.HRR131895,fig.HRR131897,fig.HRR131901,fig.HRR131904,fig.HRR131906,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt2.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS4181123,fig.SRS4181126,fig.SRS7727466,fig.SRS7727467,fig.SRS12015709,fig.SRS3822680,fig.SRS3822682,fig.SRS3822683,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt3.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS3822686,fig.LZ011,fig.SRS12015710,fig.LZ009,fig.SRS12015708,fig.SRS5086057,fig.SRS5086058,fig.LZ005,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt4.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.LZ008,fig.SRS5086059,fig.SRS5086060,fig.SRS5086061,fig.SRS5086062,fig.SRS5086063,fig.SRS5086064,fig.SRS3065428,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt5.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.LZ016,fig.SRS9921715,fig.SRS5883824,fig.SRS5883825,fig.SRS5883826,fig.SRS9921718,fig.SRS9921716,fig.SRS9921717,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt6.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS3065429,fig.SRS3065430,fig.SRS5787712,fig.SRS5883812,fig.SRS5883813,fig.SRS5883814,fig.LZ003,fig.LZ007,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt7.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS6959440,fig.SRS6959441,fig.SRS2823407,fig.SRS2823409,fig.OES026004,fig.SRS5883810,fig.SRS5883811,fig.SRS4181127,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt8.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS4181130,fig.SRS2823408,fig.SRS5086065,fig.SRS5086066,fig.SRS6959442,fig.SRS9921722,fig.SRS9921724,fig.SRS9921719,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_human_time_course_pt9.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS9921720,fig.SRS9921721,fig.SRS9921723,fig.SRS9921726,fig.SRS9921725,ncol=3,nrow=2,common.legend=TRUE,legend='bottom')
dev.off()