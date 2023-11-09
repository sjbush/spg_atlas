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

# load query ERS3000379 (5 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS3000379<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS3000379/R_SCT/ERS3000379.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS3000379) %in% mouse_genes$Symbol
ERS3000379_v2 <- ERS3000379[tmp_idx,]
num<-length(rownames(ERS3000379_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS3000379_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS3000379_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS3000379_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS3000379_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS3000379_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS3000379<-ERS3000379_v2
DefaultAssay(ERS3000379) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS3000379, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS3000379 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS3000379, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS3000379$predicted.cluster <- factor(ERS3000379$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS3000379 <- DimPlot(ERS3000379, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 days\n(ERS3000379)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS3000379$predicted.cluster))
predicted_sscs.ERS3000379 <- subset(ERS3000379@meta.data,((ERS3000379@meta.data$predicted.cluster != 3) & (ERS3000379@meta.data$predicted.cluster != 7) & (ERS3000379@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS3000379,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS3000379.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS3000379<-NULL

# load query ERS3000380 (5 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS3000380<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS3000380/R_SCT/ERS3000380.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS3000380) %in% mouse_genes$Symbol
ERS3000380_v2 <- ERS3000380[tmp_idx,]
num<-length(rownames(ERS3000380_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS3000380_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS3000380_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS3000380_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS3000380_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS3000380_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS3000380<-ERS3000380_v2
DefaultAssay(ERS3000380) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS3000380, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS3000380 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS3000380, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS3000380$predicted.cluster <- factor(ERS3000380$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS3000380 <- DimPlot(ERS3000380, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 days\n(ERS3000380)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS3000380$predicted.cluster))
predicted_sscs.ERS3000380 <- subset(ERS3000380@meta.data,((ERS3000380@meta.data$predicted.cluster != 3) & (ERS3000380@meta.data$predicted.cluster != 7) & (ERS3000380@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS3000380,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS3000380.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS3000380<-NULL

# load query SRS3990943 (6 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990943<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990943/R_SCT/SRS3990943.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990943) %in% mouse_genes$Symbol
SRS3990943_v2 <- SRS3990943[tmp_idx,]
num<-length(rownames(SRS3990943_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990943_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990943_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990943_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990943_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990943_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990943<-SRS3990943_v2
DefaultAssay(SRS3990943) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990943, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990943 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990943, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990943$predicted.cluster <- factor(SRS3990943$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990943 <- DimPlot(SRS3990943, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("6 days\n(SRS3990943)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990943$predicted.cluster))
predicted_sscs.SRS3990943 <- subset(SRS3990943@meta.data,((SRS3990943@meta.data$predicted.cluster != 3) & (SRS3990943@meta.data$predicted.cluster != 7) & (SRS3990943@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990943,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990943.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990943<-NULL

# load query ERS2575682 (10 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575682<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575682/R_SCT/ERS2575682.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575682) %in% mouse_genes$Symbol
ERS2575682_v2 <- ERS2575682[tmp_idx,]
num<-length(rownames(ERS2575682_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575682_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575682_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575682_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575682_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575682_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575682<-ERS2575682_v2
DefaultAssay(ERS2575682) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575682, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575682 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575682, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575682$predicted.cluster <- factor(ERS2575682$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575682 <- DimPlot(ERS2575682, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("10 days\n(ERS2575682)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575682$predicted.cluster))
predicted_sscs.ERS2575682 <- subset(ERS2575682@meta.data,((ERS2575682@meta.data$predicted.cluster != 3) & (ERS2575682@meta.data$predicted.cluster != 7) & (ERS2575682@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575682,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575682.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575682<-NULL

# load query SRS3990942 (14 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990942<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990942/R_SCT/SRS3990942.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990942) %in% mouse_genes$Symbol
SRS3990942_v2 <- SRS3990942[tmp_idx,]
num<-length(rownames(SRS3990942_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990942_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990942_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990942_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990942_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990942_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990942<-SRS3990942_v2
DefaultAssay(SRS3990942) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990942, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990942 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990942, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990942$predicted.cluster <- factor(SRS3990942$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990942 <- DimPlot(SRS3990942, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("14 days\n(SRS3990942)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990942$predicted.cluster))
predicted_sscs.SRS3990942 <- subset(SRS3990942@meta.data,((SRS3990942@meta.data$predicted.cluster != 3) & (SRS3990942@meta.data$predicted.cluster != 7) & (SRS3990942@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990942,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990942.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990942<-NULL

# load query ERS2575686 (15 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575686<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575686/R_SCT/ERS2575686.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575686) %in% mouse_genes$Symbol
ERS2575686_v2 <- ERS2575686[tmp_idx,]
num<-length(rownames(ERS2575686_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575686_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575686_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575686_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575686_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575686_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575686<-ERS2575686_v2
DefaultAssay(ERS2575686) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575686, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575686 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575686, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575686$predicted.cluster <- factor(ERS2575686$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575686 <- DimPlot(ERS2575686, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("15 days\n(ERS2575686)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575686$predicted.cluster))
predicted_sscs.ERS2575686 <- subset(ERS2575686@meta.data,((ERS2575686@meta.data$predicted.cluster != 3) & (ERS2575686@meta.data$predicted.cluster != 7) & (ERS2575686@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575686,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575686.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575686<-NULL

# load query SRS3990944 (18 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990944<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990944/R_SCT/SRS3990944.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990944) %in% mouse_genes$Symbol
SRS3990944_v2 <- SRS3990944[tmp_idx,]
num<-length(rownames(SRS3990944_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990944_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990944_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990944_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990944_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990944_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990944<-SRS3990944_v2
DefaultAssay(SRS3990944) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990944, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990944 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990944, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990944$predicted.cluster <- factor(SRS3990944$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990944 <- DimPlot(SRS3990944, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("18 days\n(SRS3990944)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990944$predicted.cluster))
predicted_sscs.SRS3990944 <- subset(SRS3990944@meta.data,((SRS3990944@meta.data$predicted.cluster != 3) & (SRS3990944@meta.data$predicted.cluster != 7) & (SRS3990944@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990944,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990944.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990944<-NULL

# load query SRS3990945 (18 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990945<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990945/R_SCT/SRS3990945.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990945) %in% mouse_genes$Symbol
SRS3990945_v2 <- SRS3990945[tmp_idx,]
num<-length(rownames(SRS3990945_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990945_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990945_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990945_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990945_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990945_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990945<-SRS3990945_v2
DefaultAssay(SRS3990945) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990945, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990945 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990945, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990945$predicted.cluster <- factor(SRS3990945$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990945 <- DimPlot(SRS3990945, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("18 days\n(SRS3990945)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990945$predicted.cluster))
predicted_sscs.SRS3990945 <- subset(SRS3990945@meta.data,((SRS3990945@meta.data$predicted.cluster != 3) & (SRS3990945@meta.data$predicted.cluster != 7) & (SRS3990945@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990945,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990945.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990945<-NULL

# load query ERS2575683 (20 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575683<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575683/R_SCT/ERS2575683.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575683) %in% mouse_genes$Symbol
ERS2575683_v2 <- ERS2575683[tmp_idx,]
num<-length(rownames(ERS2575683_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575683_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575683_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575683_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575683_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575683_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575683<-ERS2575683_v2
DefaultAssay(ERS2575683) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575683, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575683 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575683, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575683$predicted.cluster <- factor(ERS2575683$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575683 <- DimPlot(ERS2575683, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("20 days\n(ERS2575683)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575683$predicted.cluster))
predicted_sscs.ERS2575683 <- subset(ERS2575683@meta.data,((ERS2575683@meta.data$predicted.cluster != 3) & (ERS2575683@meta.data$predicted.cluster != 7) & (ERS2575683@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575683,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575683.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575683<-NULL

# load query SRS3990946 (25 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990946<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990946/R_SCT/SRS3990946.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990946) %in% mouse_genes$Symbol
SRS3990946_v2 <- SRS3990946[tmp_idx,]
num<-length(rownames(SRS3990946_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990946_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990946_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990946_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990946_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990946_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990946<-SRS3990946_v2
DefaultAssay(SRS3990946) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990946, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990946 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990946, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990946$predicted.cluster <- factor(SRS3990946$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990946 <- DimPlot(SRS3990946, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("25 days\n(SRS3990946)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990946$predicted.cluster))
predicted_sscs.SRS3990946 <- subset(SRS3990946@meta.data,((SRS3990946@meta.data$predicted.cluster != 3) & (SRS3990946@meta.data$predicted.cluster != 7) & (SRS3990946@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990946,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990946.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990946<-NULL

# load query ERS2575687 (25 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575687<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575687/R_SCT/ERS2575687.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575687) %in% mouse_genes$Symbol
ERS2575687_v2 <- ERS2575687[tmp_idx,]
num<-length(rownames(ERS2575687_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575687_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575687_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575687_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575687_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575687_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575687<-ERS2575687_v2
DefaultAssay(ERS2575687) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575687, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575687 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575687, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575687$predicted.cluster <- factor(ERS2575687$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575687 <- DimPlot(ERS2575687, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("25 days\n(ERS2575687)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575687$predicted.cluster))
predicted_sscs.ERS2575687 <- subset(ERS2575687@meta.data,((ERS2575687@meta.data$predicted.cluster != 3) & (ERS2575687@meta.data$predicted.cluster != 7) & (ERS2575687@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575687,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575687.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575687<-NULL

# load query SRS3990947 (30 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990947<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990947/R_SCT/SRS3990947.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990947) %in% mouse_genes$Symbol
SRS3990947_v2 <- SRS3990947[tmp_idx,]
num<-length(rownames(SRS3990947_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990947_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990947_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990947_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990947_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990947_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990947<-SRS3990947_v2
DefaultAssay(SRS3990947) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990947, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990947 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990947, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990947$predicted.cluster <- factor(SRS3990947$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990947 <- DimPlot(SRS3990947, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("30 days\n(SRS3990947)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990947$predicted.cluster))
predicted_sscs.SRS3990947 <- subset(SRS3990947@meta.data,((SRS3990947@meta.data$predicted.cluster != 3) & (SRS3990947@meta.data$predicted.cluster != 7) & (SRS3990947@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990947,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990947.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990947<-NULL

# load query ERS2575684 (30 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575684<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575684/R_SCT/ERS2575684.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575684) %in% mouse_genes$Symbol
ERS2575684_v2 <- ERS2575684[tmp_idx,]
num<-length(rownames(ERS2575684_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575684_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575684_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575684_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575684_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575684_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575684<-ERS2575684_v2
DefaultAssay(ERS2575684) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575684, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575684 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575684, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575684$predicted.cluster <- factor(ERS2575684$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575684 <- DimPlot(ERS2575684, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("30 days\n(ERS2575684)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575684$predicted.cluster))
predicted_sscs.ERS2575684 <- subset(ERS2575684@meta.data,((ERS2575684@meta.data$predicted.cluster != 3) & (ERS2575684@meta.data$predicted.cluster != 7) & (ERS2575684@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575684,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575684.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575684<-NULL

# load query ERS2575685 (35 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575685<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575685/R_SCT/ERS2575685.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575685) %in% mouse_genes$Symbol
ERS2575685_v2 <- ERS2575685[tmp_idx,]
num<-length(rownames(ERS2575685_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575685_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575685_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575685_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575685_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575685_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575685<-ERS2575685_v2
DefaultAssay(ERS2575685) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575685, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575685 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575685, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575685$predicted.cluster <- factor(ERS2575685$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575685 <- DimPlot(ERS2575685, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("35 days\n(ERS2575685)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575685$predicted.cluster))
predicted_sscs.ERS2575685 <- subset(ERS2575685@meta.data,((ERS2575685@meta.data$predicted.cluster != 3) & (ERS2575685@meta.data$predicted.cluster != 7) & (ERS2575685@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575685,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575685.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575685<-NULL

# load query SRS3990948 (56 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990948<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990948/R_SCT/SRS3990948.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990948) %in% mouse_genes$Symbol
SRS3990948_v2 <- SRS3990948[tmp_idx,]
num<-length(rownames(SRS3990948_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990948_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990948_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990948_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990948_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990948_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990948<-SRS3990948_v2
DefaultAssay(SRS3990948) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990948, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990948 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990948, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990948$predicted.cluster <- factor(SRS3990948$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990948 <- DimPlot(SRS3990948, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("56 days\n(SRS3990948)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990948$predicted.cluster))
predicted_sscs.SRS3990948 <- subset(SRS3990948@meta.data,((SRS3990948@meta.data$predicted.cluster != 3) & (SRS3990948@meta.data$predicted.cluster != 7) & (SRS3990948@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990948,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990948.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990948<-NULL

# load query SRS3990949 (56 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3990949<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3990949/R_SCT/SRS3990949.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3990949) %in% mouse_genes$Symbol
SRS3990949_v2 <- SRS3990949[tmp_idx,]
num<-length(rownames(SRS3990949_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3990949_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3990949_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3990949_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3990949_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3990949_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3990949<-SRS3990949_v2
DefaultAssay(SRS3990949) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3990949, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3990949 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3990949, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3990949$predicted.cluster <- factor(SRS3990949$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3990949 <- DimPlot(SRS3990949, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("56 days\n(SRS3990949)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3990949$predicted.cluster))
predicted_sscs.SRS3990949 <- subset(SRS3990949@meta.data,((SRS3990949@meta.data$predicted.cluster != 3) & (SRS3990949@meta.data$predicted.cluster != 7) & (SRS3990949@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3990949,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3990949.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3990949<-NULL

# load query ERS2575688 (61 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575688<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575688/R_SCT/ERS2575688.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575688) %in% mouse_genes$Symbol
ERS2575688_v2 <- ERS2575688[tmp_idx,]
num<-length(rownames(ERS2575688_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575688_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575688_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575688_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575688_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575688_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575688<-ERS2575688_v2
DefaultAssay(ERS2575688) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575688, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575688 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575688, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575688$predicted.cluster <- factor(ERS2575688$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575688 <- DimPlot(ERS2575688, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("61 days\n(ERS2575688)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575688$predicted.cluster))
predicted_sscs.ERS2575688 <- subset(ERS2575688@meta.data,((ERS2575688@meta.data$predicted.cluster != 3) & (ERS2575688@meta.data$predicted.cluster != 7) & (ERS2575688@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575688,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575688.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575688<-NULL

# load query ERS2575689 (61 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575689<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575689/R_SCT/ERS2575689.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575689) %in% mouse_genes$Symbol
ERS2575689_v2 <- ERS2575689[tmp_idx,]
num<-length(rownames(ERS2575689_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575689_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575689_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575689_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575689_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575689_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575689<-ERS2575689_v2
DefaultAssay(ERS2575689) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575689, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575689 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575689, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575689$predicted.cluster <- factor(ERS2575689$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575689 <- DimPlot(ERS2575689, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("61 days\n(ERS2575689)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575689$predicted.cluster))
predicted_sscs.ERS2575689 <- subset(ERS2575689@meta.data,((ERS2575689@meta.data$predicted.cluster != 3) & (ERS2575689@meta.data$predicted.cluster != 7) & (ERS2575689@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575689,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575689.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575689<-NULL

# load query ERS2575690 (62 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575690<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575690/R_SCT/ERS2575690.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575690) %in% mouse_genes$Symbol
ERS2575690_v2 <- ERS2575690[tmp_idx,]
num<-length(rownames(ERS2575690_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575690_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575690_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575690_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575690_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575690_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575690<-ERS2575690_v2
DefaultAssay(ERS2575690) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575690, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575690 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575690, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575690$predicted.cluster <- factor(ERS2575690$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575690 <- DimPlot(ERS2575690, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("62 days\n(ERS2575690)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575690$predicted.cluster))
predicted_sscs.ERS2575690 <- subset(ERS2575690@meta.data,((ERS2575690@meta.data$predicted.cluster != 3) & (ERS2575690@meta.data$predicted.cluster != 7) & (ERS2575690@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575690,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575690.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575690<-NULL

# load query ERS2575678 (63 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575678<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575678/R_SCT/ERS2575678.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575678) %in% mouse_genes$Symbol
ERS2575678_v2 <- ERS2575678[tmp_idx,]
num<-length(rownames(ERS2575678_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575678_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575678_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575678_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575678_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575678_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575678<-ERS2575678_v2
DefaultAssay(ERS2575678) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575678, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575678 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575678, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575678$predicted.cluster <- factor(ERS2575678$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575678 <- DimPlot(ERS2575678, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("63 days\n(ERS2575678)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575678$predicted.cluster))
predicted_sscs.ERS2575678 <- subset(ERS2575678@meta.data,((ERS2575678@meta.data$predicted.cluster != 3) & (ERS2575678@meta.data$predicted.cluster != 7) & (ERS2575678@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575678,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575678.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575678<-NULL

# load query ERS2575679 (63 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575679<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575679/R_SCT/ERS2575679.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575679) %in% mouse_genes$Symbol
ERS2575679_v2 <- ERS2575679[tmp_idx,]
num<-length(rownames(ERS2575679_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575679_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575679_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575679_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575679_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575679_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575679<-ERS2575679_v2
DefaultAssay(ERS2575679) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575679, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575679 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575679, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575679$predicted.cluster <- factor(ERS2575679$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575679 <- DimPlot(ERS2575679, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("63 days\n(ERS2575679)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575679$predicted.cluster))
predicted_sscs.ERS2575679 <- subset(ERS2575679@meta.data,((ERS2575679@meta.data$predicted.cluster != 3) & (ERS2575679@meta.data$predicted.cluster != 7) & (ERS2575679@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575679,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575679.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575679<-NULL

# load query SRS3097511 (63 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3097511<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3097511/R_SCT/SRS3097511.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3097511) %in% mouse_genes$Symbol
SRS3097511_v2 <- SRS3097511[tmp_idx,]
num<-length(rownames(SRS3097511_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3097511_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3097511_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3097511_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3097511_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3097511_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3097511<-SRS3097511_v2
DefaultAssay(SRS3097511) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3097511, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3097511 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3097511, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3097511$predicted.cluster <- factor(SRS3097511$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3097511 <- DimPlot(SRS3097511, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("63 days\n(SRS3097511)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3097511$predicted.cluster))
predicted_sscs.SRS3097511 <- subset(SRS3097511@meta.data,((SRS3097511@meta.data$predicted.cluster != 3) & (SRS3097511@meta.data$predicted.cluster != 7) & (SRS3097511@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3097511,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3097511.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3097511<-NULL

# load query SRS3097512 (63 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3097512<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3097512/R_SCT/SRS3097512.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3097512) %in% mouse_genes$Symbol
SRS3097512_v2 <- SRS3097512[tmp_idx,]
num<-length(rownames(SRS3097512_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3097512_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3097512_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3097512_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3097512_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3097512_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3097512<-SRS3097512_v2
DefaultAssay(SRS3097512) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3097512, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3097512 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3097512, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3097512$predicted.cluster <- factor(SRS3097512$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3097512 <- DimPlot(SRS3097512, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("63 days\n(SRS3097512)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3097512$predicted.cluster))
predicted_sscs.SRS3097512 <- subset(SRS3097512@meta.data,((SRS3097512@meta.data$predicted.cluster != 3) & (SRS3097512@meta.data$predicted.cluster != 7) & (SRS3097512@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3097512,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3097512.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3097512<-NULL

# load query SRS3097514 (63 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3097514<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3097514/R_SCT/SRS3097514.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3097514) %in% mouse_genes$Symbol
SRS3097514_v2 <- SRS3097514[tmp_idx,]
num<-length(rownames(SRS3097514_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3097514_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3097514_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3097514_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3097514_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3097514_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3097514<-SRS3097514_v2
DefaultAssay(SRS3097514) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3097514, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3097514 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3097514, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3097514$predicted.cluster <- factor(SRS3097514$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3097514 <- DimPlot(SRS3097514, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("63 days\n(SRS3097514)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3097514$predicted.cluster))
predicted_sscs.SRS3097514 <- subset(SRS3097514@meta.data,((SRS3097514@meta.data$predicted.cluster != 3) & (SRS3097514@meta.data$predicted.cluster != 7) & (SRS3097514@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3097514,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3097514.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3097514<-NULL

# load query SRS3097515 (63 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3097515<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3097515/R_SCT/SRS3097515.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3097515) %in% mouse_genes$Symbol
SRS3097515_v2 <- SRS3097515[tmp_idx,]
num<-length(rownames(SRS3097515_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3097515_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3097515_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3097515_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3097515_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3097515_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3097515<-SRS3097515_v2
DefaultAssay(SRS3097515) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3097515, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3097515 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3097515, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3097515$predicted.cluster <- factor(SRS3097515$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3097515 <- DimPlot(SRS3097515, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("63 days\n(SRS3097515)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3097515$predicted.cluster))
predicted_sscs.SRS3097515 <- subset(SRS3097515@meta.data,((SRS3097515@meta.data$predicted.cluster != 3) & (SRS3097515@meta.data$predicted.cluster != 7) & (SRS3097515@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3097515,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3097515.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3097515<-NULL

# load query ERS2575676 (67 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575676<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575676/R_SCT/ERS2575676.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575676) %in% mouse_genes$Symbol
ERS2575676_v2 <- ERS2575676[tmp_idx,]
num<-length(rownames(ERS2575676_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575676_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575676_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575676_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575676_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575676_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575676<-ERS2575676_v2
DefaultAssay(ERS2575676) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575676, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575676 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575676, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575676$predicted.cluster <- factor(ERS2575676$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575676 <- DimPlot(ERS2575676, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("67 days\n(ERS2575676)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575676$predicted.cluster))
predicted_sscs.ERS2575676 <- subset(ERS2575676@meta.data,((ERS2575676@meta.data$predicted.cluster != 3) & (ERS2575676@meta.data$predicted.cluster != 7) & (ERS2575676@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575676,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575676.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575676<-NULL

# load query ERS2575677 (67 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575677<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575677/R_SCT/ERS2575677.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575677) %in% mouse_genes$Symbol
ERS2575677_v2 <- ERS2575677[tmp_idx,]
num<-length(rownames(ERS2575677_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575677_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575677_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575677_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575677_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575677_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575677<-ERS2575677_v2
DefaultAssay(ERS2575677) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575677, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575677 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575677, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575677$predicted.cluster <- factor(ERS2575677$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575677 <- DimPlot(ERS2575677, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("67 days\n(ERS2575677)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575677$predicted.cluster))
predicted_sscs.ERS2575677 <- subset(ERS2575677@meta.data,((ERS2575677@meta.data$predicted.cluster != 3) & (ERS2575677@meta.data$predicted.cluster != 7) & (ERS2575677@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575677,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575677.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575677<-NULL

# load query ERS2575680 (67 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575680<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575680/R_SCT/ERS2575680.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575680) %in% mouse_genes$Symbol
ERS2575680_v2 <- ERS2575680[tmp_idx,]
num<-length(rownames(ERS2575680_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575680_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575680_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575680_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575680_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575680_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575680<-ERS2575680_v2
DefaultAssay(ERS2575680) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575680, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575680 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575680, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575680$predicted.cluster <- factor(ERS2575680$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575680 <- DimPlot(ERS2575680, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("67 days\n(ERS2575680)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575680$predicted.cluster))
predicted_sscs.ERS2575680 <- subset(ERS2575680@meta.data,((ERS2575680@meta.data$predicted.cluster != 3) & (ERS2575680@meta.data$predicted.cluster != 7) & (ERS2575680@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575680,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575680.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575680<-NULL

# load query ERS2575681 (67 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
ERS2575681<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/ERS2575681/R_SCT/ERS2575681.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(ERS2575681) %in% mouse_genes$Symbol
ERS2575681_v2 <- ERS2575681[tmp_idx,]
num<-length(rownames(ERS2575681_v2))
genetable_1 <- data.frame(mSymbol=rownames(ERS2575681_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(ERS2575681_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(ERS2575681_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(ERS2575681_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(ERS2575681_v2@assays$SCT@data) <- genetable_2$HumanSymbol
ERS2575681<-ERS2575681_v2
DefaultAssay(ERS2575681) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = ERS2575681, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
ERS2575681 <- MapQuery(anchorset = adult.anchors, reference = adult, query = ERS2575681, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
ERS2575681$predicted.cluster <- factor(ERS2575681$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.ERS2575681 <- DimPlot(ERS2575681, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("67 days\n(ERS2575681)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(ERS2575681$predicted.cluster))
predicted_sscs.ERS2575681 <- subset(ERS2575681@meta.data,((ERS2575681@meta.data$predicted.cluster != 3) & (ERS2575681@meta.data$predicted.cluster != 7) & (ERS2575681@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.ERS2575681,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.ERS2575681.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
ERS2575681<-NULL

# load query SRS3189008 (3 months 1 day), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189008<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189008/R_SCT/SRS3189008.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189008) %in% mouse_genes$Symbol
SRS3189008_v2 <- SRS3189008[tmp_idx,]
num<-length(rownames(SRS3189008_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189008_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189008_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189008_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189008_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189008_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189008<-SRS3189008_v2
DefaultAssay(SRS3189008) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189008, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189008 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189008, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189008$predicted.cluster <- factor(SRS3189008$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189008 <- DimPlot(SRS3189008, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("3 months 1 day\n(SRS3189008)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189008$predicted.cluster))
predicted_sscs.SRS3189008 <- subset(SRS3189008@meta.data,((SRS3189008@meta.data$predicted.cluster != 3) & (SRS3189008@meta.data$predicted.cluster != 7) & (SRS3189008@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189008,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189008.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189008<-NULL

# load query SRS3189011 (3 months 2 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189011<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189011/R_SCT/SRS3189011.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189011) %in% mouse_genes$Symbol
SRS3189011_v2 <- SRS3189011[tmp_idx,]
num<-length(rownames(SRS3189011_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189011_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189011_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189011_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189011_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189011_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189011<-SRS3189011_v2
DefaultAssay(SRS3189011) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189011, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189011 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189011, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189011$predicted.cluster <- factor(SRS3189011$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189011 <- DimPlot(SRS3189011, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("3 months 2 days\n(SRS3189011)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189011$predicted.cluster))
predicted_sscs.SRS3189011 <- subset(SRS3189011@meta.data,((SRS3189011@meta.data$predicted.cluster != 3) & (SRS3189011@meta.data$predicted.cluster != 7) & (SRS3189011@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189011,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189011.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189011<-NULL

# load query SRS3189004 (3 month 21 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189004<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189004/R_SCT/SRS3189004.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189004) %in% mouse_genes$Symbol
SRS3189004_v2 <- SRS3189004[tmp_idx,]
num<-length(rownames(SRS3189004_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189004_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189004_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189004_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189004_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189004_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189004<-SRS3189004_v2
DefaultAssay(SRS3189004) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189004, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189004 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189004, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189004$predicted.cluster <- factor(SRS3189004$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189004 <- DimPlot(SRS3189004, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("3 months 21 days\n(SRS3189004)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189004$predicted.cluster))
predicted_sscs.SRS3189004 <- subset(SRS3189004@meta.data,((SRS3189004@meta.data$predicted.cluster != 3) & (SRS3189004@meta.data$predicted.cluster != 7) & (SRS3189004@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189004,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189004.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189004<-NULL

# load query SRS3189007 (3 month 27 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189007<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189007/R_SCT/SRS3189007.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189007) %in% mouse_genes$Symbol
SRS3189007_v2 <- SRS3189007[tmp_idx,]
num<-length(rownames(SRS3189007_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189007_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189007_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189007_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189007_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189007_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189007<-SRS3189007_v2
DefaultAssay(SRS3189007) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189007, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189007 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189007, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189007$predicted.cluster <- factor(SRS3189007$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189007 <- DimPlot(SRS3189007, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("3 months 27 days\n(SRS3189007)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189007$predicted.cluster))
predicted_sscs.SRS3189007 <- subset(SRS3189007@meta.data,((SRS3189007@meta.data$predicted.cluster != 3) & (SRS3189007@meta.data$predicted.cluster != 7) & (SRS3189007@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189007,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189007.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189007<-NULL

# load query SRS3189006 (4 months), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189006<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189006/R_SCT/SRS3189006.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189006) %in% mouse_genes$Symbol
SRS3189006_v2 <- SRS3189006[tmp_idx,]
num<-length(rownames(SRS3189006_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189006_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189006_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189006_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189006_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189006_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189006<-SRS3189006_v2
DefaultAssay(SRS3189006) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189006, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189006 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189006, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189006$predicted.cluster <- factor(SRS3189006$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189006 <- DimPlot(SRS3189006, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("4 months\n(SRS3189006)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189006$predicted.cluster))
predicted_sscs.SRS3189006 <- subset(SRS3189006@meta.data,((SRS3189006@meta.data$predicted.cluster != 3) & (SRS3189006@meta.data$predicted.cluster != 7) & (SRS3189006@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189006,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189006.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189006<-NULL

# load query SRS3189009 (4 months), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189009<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189009/R_SCT/SRS3189009.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189009) %in% mouse_genes$Symbol
SRS3189009_v2 <- SRS3189009[tmp_idx,]
num<-length(rownames(SRS3189009_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189009_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189009_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189009_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189009_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189009_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189009<-SRS3189009_v2
DefaultAssay(SRS3189009) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189009, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189009 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189009, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189009$predicted.cluster <- factor(SRS3189009$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189009 <- DimPlot(SRS3189009, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("4 months\n(SRS3189009)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189009$predicted.cluster))
predicted_sscs.SRS3189009 <- subset(SRS3189009@meta.data,((SRS3189009@meta.data$predicted.cluster != 3) & (SRS3189009@meta.data$predicted.cluster != 7) & (SRS3189009@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189009,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189009.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189009<-NULL

# load query SRS3189026 (4 months 2 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189026<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189026/R_SCT/SRS3189026.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189026) %in% mouse_genes$Symbol
SRS3189026_v2 <- SRS3189026[tmp_idx,]
num<-length(rownames(SRS3189026_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189026_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189026_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189026_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189026_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189026_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189026<-SRS3189026_v2
DefaultAssay(SRS3189026) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189026, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189026 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189026, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189026$predicted.cluster <- factor(SRS3189026$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189026 <- DimPlot(SRS3189026, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("4 months 2 days\n(SRS3189026)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189026$predicted.cluster))
predicted_sscs.SRS3189026 <- subset(SRS3189026@meta.data,((SRS3189026@meta.data$predicted.cluster != 3) & (SRS3189026@meta.data$predicted.cluster != 7) & (SRS3189026@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189026,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189026.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189026<-NULL

# load query SRS3189010 (5 months 5 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189010<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189010/R_SCT/SRS3189010.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189010) %in% mouse_genes$Symbol
SRS3189010_v2 <- SRS3189010[tmp_idx,]
num<-length(rownames(SRS3189010_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189010_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189010_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189010_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189010_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189010_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189010<-SRS3189010_v2
DefaultAssay(SRS3189010) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189010, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189010 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189010, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189010$predicted.cluster <- factor(SRS3189010$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189010 <- DimPlot(SRS3189010, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 months 5 days\n(SRS3189010)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189010$predicted.cluster))
predicted_sscs.SRS3189010 <- subset(SRS3189010@meta.data,((SRS3189010@meta.data$predicted.cluster != 3) & (SRS3189010@meta.data$predicted.cluster != 7) & (SRS3189010@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189010,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189010.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189010<-NULL

# load query SRS3189005 (5 months 7 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189005<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189005/R_SCT/SRS3189005.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189005) %in% mouse_genes$Symbol
SRS3189005_v2 <- SRS3189005[tmp_idx,]
num<-length(rownames(SRS3189005_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189005_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189005_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189005_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189005_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189005_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189005<-SRS3189005_v2
DefaultAssay(SRS3189005) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189005, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189005 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189005, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189005$predicted.cluster <- factor(SRS3189005$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189005 <- DimPlot(SRS3189005, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("5 months 7 days\n(SRS3189005)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189005$predicted.cluster))
predicted_sscs.SRS3189005 <- subset(SRS3189005@meta.data,((SRS3189005@meta.data$predicted.cluster != 3) & (SRS3189005@meta.data$predicted.cluster != 7) & (SRS3189005@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189005,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189005.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189005<-NULL

# load query SRS3189013 (7 months 21 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189013<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189013/R_SCT/SRS3189013.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189013) %in% mouse_genes$Symbol
SRS3189013_v2 <- SRS3189013[tmp_idx,]
num<-length(rownames(SRS3189013_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189013_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189013_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189013_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189013_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189013_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189013<-SRS3189013_v2
DefaultAssay(SRS3189013) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189013, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189013 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189013, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189013$predicted.cluster <- factor(SRS3189013$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189013 <- DimPlot(SRS3189013, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("7 months 21 days\n(SRS3189013)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189013$predicted.cluster))
predicted_sscs.SRS3189013 <- subset(SRS3189013@meta.data,((SRS3189013@meta.data$predicted.cluster != 3) & (SRS3189013@meta.data$predicted.cluster != 7) & (SRS3189013@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189013,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189013.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189013<-NULL

# load query SRS3189012 (9 months 14 days), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3189012<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3189012/R_SCT/SRS3189012.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3189012) %in% mouse_genes$Symbol
SRS3189012_v2 <- SRS3189012[tmp_idx,]
num<-length(rownames(SRS3189012_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3189012_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3189012_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3189012_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3189012_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3189012_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3189012<-SRS3189012_v2
DefaultAssay(SRS3189012) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3189012, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3189012 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3189012, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("0","1","2","3","4","5","6","7","8","9") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3189012$predicted.cluster <- factor(SRS3189012$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3189012 <- DimPlot(SRS3189012, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("9 months 14 days\n(SRS3189012)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3189012$predicted.cluster))
predicted_sscs.SRS3189012 <- subset(SRS3189012@meta.data,((SRS3189012@meta.data$predicted.cluster != 3) & (SRS3189012@meta.data$predicted.cluster != 7) & (SRS3189012@meta.data$predicted.cluster != 9)))
write.table(predicted_sscs.SRS3189012,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS3189012.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS3189012<-NULL

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_mouse_pt1.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.ERS3000379,fig.ERS3000380,fig.SRS3990943,fig.ERS2575682,fig.SRS3990942,fig.ERS2575686,fig.SRS3990944,fig.SRS3990945,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_mouse_pt2.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.ERS2575683,fig.SRS3990946,fig.ERS2575687,fig.SRS3990947,fig.ERS2575684,fig.ERS2575685,fig.SRS3990948,fig.SRS3990949,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_mouse_pt3.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.ERS2575688,fig.ERS2575689,fig.ERS2575690,fig.ERS2575678,fig.ERS2575679,fig.SRS3097511,fig.SRS3097512,fig.SRS3097514,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_mouse_pt4.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS3097515,fig.ERS2575676,fig.ERS2575677,fig.ERS2575680,fig.ERS2575681,fig.SRS3189008,fig.SRS3189011,fig.SRS3189004,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_mouse_pt5.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS3189007,fig.SRS3189006,fig.SRS3189009,fig.SRS3189026,fig.SRS3189010,fig.SRS3189005,fig.SRS3189013,fig.SRS3189012,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()