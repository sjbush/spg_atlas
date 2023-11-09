## PROJECTIONS: of sheep, pig, buffalo and cynomolgus macaque samples onto the states 0-4 UMAP. This is a gating strategy to distinguish somatic from germline cells, before re-projecting the latter onto the states 0-4 UMAP
## IMPORTANT: use "module add R-cbrg/202210" to ensure you use v4.1.1 of Seurat
## NOTE: the X_to_human_symbol_lookup.txt files are created by 16pre2.create_SpeciesX_to_human_symbol_lookup_table.pl

library(Seurat)
library(ggpubr)
library(tidyverse)
theme_set(theme_bw())

# load reference: adult human states 0-4 data
hSSC <- readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(hSSC)<-hSSC$integrated_snn_res.1.1
new.cluster.ids <- c(
"state 0A/1",
"SPG",
"state 0",
"leptotene",
"state 0B",
"zygotene",
"diff SPG"
)
names(new.cluster.ids) <- levels(hSSC)
hSSC <- RenameIdents(hSSC, new.cluster.ids)
hSSC@meta.data$ssc_state <- Idents(hSSC)

# to re-create the human UMAP, we first need to quantify the content of the elbow plot, implementing code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- hSSC[["pca"]]@stdev / sum(hSSC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let us choose the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data. We need this parameter to input into RunUMAP
prin_comp <- min(component1, component2)

# we have previously obtained the list of predicted cluster identities after projecting each human sample onto the adult whole-testes atlas (which had been clustered at resolution 0.15)
# we had subset that list to contain only those cells which map to the germline - not somatic - clusters
# we will now project only germ cells onto the SSC atlas
# the way to do this is NOT to subset the Seurat object as that will leave you with too few cells (the k.score and k.filter parameters of FindTransferAnchors will cause errors); rather, we tell DimPlot to only print certain cells
predicted_sscs.SRS7528324<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS7528324.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9029393<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9029393.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS11257878<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS11257878.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS11257879<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS11257879.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.infant<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.infant.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.juvenile<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.juvenile.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.adult1<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.adult1.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.adult2<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.adult2.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')

cells_to_use.SRS7528324<-predicted_sscs.SRS7528324$barcode
cells_to_use.SRS9029393<-predicted_sscs.SRS9029393$barcode
cells_to_use.SRS11257878<-predicted_sscs.SRS11257878$barcode
cells_to_use.SRS11257879<-predicted_sscs.SRS11257879$barcode
cells_to_use.infant<-predicted_sscs.infant$barcode
cells_to_use.juvenile<-predicted_sscs.juvenile$barcode
cells_to_use.adult1<-predicted_sscs.adult1$barcode
cells_to_use.adult2<-predicted_sscs.adult2$barcode

# load query SRS7528324 (sheep, 1.5 years old, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS7528324<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Ovis_aries/SRS7528324/R_SCT/SRS7528324.rds')
sheep_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/sheep_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS7528324) %in% sheep_genes$Symbol
SRS7528324_v2 <- SRS7528324[tmp_idx,]
num<-length(rownames(SRS7528324_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS7528324_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=sheep_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS7528324_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS7528324_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS7528324_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS7528324_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS7528324<-SRS7528324_v2
DefaultAssay(SRS7528324) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS7528324, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS7528324 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS7528324, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS7528324$predicted.cluster <- factor(SRS7528324$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(hSSC, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(hSSC$ssc_state), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS7528324 <- DimPlot(SRS7528324, cells = cells_to_use.SRS7528324, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("sheep\n1.5 years old") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS7528324$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS7528324@meta.data,file='states0to4.onto_which_projected_sheep_germline.SRS7528324.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS7528324<-NULL

# load query SRS9029393 (pig, 150 days old, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9029393<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Sus_scrofa/SRS9029393/R_SCT/SRS9029393.rds')
pig_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/pig_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9029393) %in% pig_genes$Symbol
SRS9029393_v2 <- SRS9029393[tmp_idx,]
num<-length(rownames(SRS9029393_v2))
genetable_1 <- data.frame(pSymbol=rownames(SRS9029393_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=pig_genes, by.x="pSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9029393_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9029393_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9029393_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9029393_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9029393<-SRS9029393_v2
DefaultAssay(SRS9029393) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9029393, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9029393 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9029393, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9029393$predicted.cluster <- factor(SRS9029393$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9029393 <- DimPlot(SRS9029393, cells = cells_to_use.SRS9029393, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("pig\n150 days old") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9029393$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9029393@meta.data,file='states0to4.onto_which_projected_pig_germline.SRS9029393.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9029393<-NULL

# load query SRS11257878 (buffalo, 3 months old, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS11257878<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Bubalus_bubalis/SRS11257878/R_SCT/SRS11257878.rds')
DefaultAssay(SRS11257878) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS11257878, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS11257878 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS11257878, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS11257878$predicted.cluster <- factor(SRS11257878$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS11257878 <- DimPlot(SRS11257878, cells = cells_to_use.SRS11257878, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("buffalo\n3 months old") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS11257878$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS11257878@meta.data,file='states0to4.onto_which_projected_buffalo_germline.SRS11257878.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS11257878<-NULL

# load query SRS11257879 (buffalo, 2 years old, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS11257879<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Bubalus_bubalis/SRS11257879/R_SCT/SRS11257879.rds')
DefaultAssay(SRS11257879) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS11257879, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS11257879 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS11257879, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS11257879$predicted.cluster <- factor(SRS11257879$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS11257879 <- DimPlot(SRS11257879, cells = cells_to_use.SRS11257879, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("buffalo\n2 years old") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS11257879$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS11257879@meta.data,file='states0to4.onto_which_projected_buffalo_germline.SRS11257879.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS11257879<-NULL

# load query infant (cynomologus macaque, 1 year), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
infant<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Macaca_fascicularis/infant/R_SCT/infant.rds')
cynomolgus_macaque_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/cynomolgus_macaque_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(infant) %in% cynomolgus_macaque_genes$Symbol
infant_v2 <- infant[tmp_idx,]
num<-length(rownames(infant_v2))
genetable_1 <- data.frame(cmSymbol=rownames(infant_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=cynomolgus_macaque_genes, by.x="cmSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(infant_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(infant_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(infant_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(infant_v2@assays$SCT@data) <- genetable_2$HumanSymbol
infant<-infant_v2
DefaultAssay(infant) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = infant, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
infant <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = infant, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
infant$predicted.cluster <- factor(infant$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.infant <- DimPlot(infant, cells = cells_to_use.infant, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("cynomolgus macaque\n1 year old") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(infant$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(infant@meta.data,file='states0to4.onto_which_projected_cynomolgus_macaque_germline.infant.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
infant<-NULL

# load query juvenile (cynomologus macaque, 2 years), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
juvenile<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Macaca_fascicularis/juvenile/R_SCT/juvenile.rds')
cynomolgus_macaque_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/cynomolgus_macaque_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(juvenile) %in% cynomolgus_macaque_genes$Symbol
juvenile_v2 <- juvenile[tmp_idx,]
num<-length(rownames(juvenile_v2))
genetable_1 <- data.frame(cmSymbol=rownames(juvenile_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=cynomolgus_macaque_genes, by.x="cmSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(juvenile_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(juvenile_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(juvenile_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(juvenile_v2@assays$SCT@data) <- genetable_2$HumanSymbol
juvenile<-juvenile_v2
DefaultAssay(juvenile) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = juvenile, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
juvenile <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = juvenile, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
juvenile$predicted.cluster <- factor(juvenile$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.juvenile <- DimPlot(juvenile, cells = cells_to_use.juvenile, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("cynomolgus macaque\n2 years old") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(juvenile$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(juvenile@meta.data,file='states0to4.onto_which_projected_cynomolgus_macaque_germline.juvenile.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
juvenile<-NULL

# load query adult1 (cynomologus macaque, > 4 years), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
adult1<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Macaca_fascicularis/adult1/R_SCT/adult1.rds')
cynomolgus_macaque_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/cynomolgus_macaque_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(adult1) %in% cynomolgus_macaque_genes$Symbol
adult1_v2 <- adult1[tmp_idx,]
num<-length(rownames(adult1_v2))
genetable_1 <- data.frame(cmSymbol=rownames(adult1_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=cynomolgus_macaque_genes, by.x="cmSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(adult1_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(adult1_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(adult1_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(adult1_v2@assays$SCT@data) <- genetable_2$HumanSymbol
adult1<-adult1_v2
DefaultAssay(adult1) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = adult1, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
adult1 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = adult1, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
adult1$predicted.cluster <- factor(adult1$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.adult1 <- DimPlot(adult1, cells = cells_to_use.adult1, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("cynomolgus macaque\n> 4 years old (Lau 2020 'adult 1')") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(adult1$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(adult1@meta.data,file='states0to4.onto_which_projected_cynomolgus_macaque_germline.adult1.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
adult1<-NULL

# load query adult2 (cynomologus macaque, > 4 years), edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
adult2<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Macaca_fascicularis/adult2/R_SCT/adult2.rds')
cynomolgus_macaque_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/cynomolgus_macaque_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(adult2) %in% cynomolgus_macaque_genes$Symbol
adult2_v2 <- adult2[tmp_idx,]
num<-length(rownames(adult2_v2))
genetable_1 <- data.frame(cmSymbol=rownames(adult2_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=cynomolgus_macaque_genes, by.x="cmSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(adult2_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(adult2_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(adult2_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(adult2_v2@assays$SCT@data) <- genetable_2$HumanSymbol
adult2<-adult2_v2
DefaultAssay(adult2) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = adult2, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
adult2 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = adult2, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
adult2$predicted.cluster <- factor(adult2$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.adult2 <- DimPlot(adult2, cells = cells_to_use.adult2, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("cynomolgus macaque\n> 4 years old (Lau 2020 'adult 2')") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(adult2$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(adult2@meta.data,file='states0to4.onto_which_projected_cynomolgus_macaque_germline.adult2.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
adult2<-NULL

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.onto_which_projected_sheep_pig_buffalo_and_cynomologus_macaque_germline.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS7528324,fig.SRS9029393,fig.SRS11257878,fig.SRS11257879,fig.infant,fig.juvenile,fig.adult1,fig.adult2,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()