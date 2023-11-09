## PROJECTIONS: of rhesus macaque samples onto the states 0-4 UMAP. This is a gating strategy to distinguish somatic from germline cells, before re-projecting the latter onto the states 0-4 UMAP
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
predicted_sscs.SRS9839846 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839846.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839847 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839847.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839848 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839848.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839849 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839849.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839850 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839850.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839851 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839851.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839852 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839852.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS9839853 <- read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected.SRS9839853.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')

cells_to_use.SRS9839846  <- predicted_sscs.SRS9839846$barcode
cells_to_use.SRS9839847 <- predicted_sscs.SRS9839847$barcode
cells_to_use.SRS9839848 <- predicted_sscs.SRS9839848$barcode
cells_to_use.SRS9839849 <- predicted_sscs.SRS9839849$barcode
cells_to_use.SRS9839850 <- predicted_sscs.SRS9839850$barcode
cells_to_use.SRS9839851 <- predicted_sscs.SRS9839851$barcode
cells_to_use.SRS9839852 <- predicted_sscs.SRS9839852$barcode
cells_to_use.SRS9839853 <- predicted_sscs.SRS9839853$barcode

# load query SRS9839846, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839846<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839846/R_SCT/SRS9839846.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839846) %in% rat_genes$Symbol
SRS9839846_v2 <- SRS9839846[tmp_idx,]
num<-length(rownames(SRS9839846_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839846_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839846_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839846_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839846_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839846_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839846<-SRS9839846_v2
DefaultAssay(SRS9839846) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839846, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839846 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839846, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839846$predicted.cluster <- factor(SRS9839846$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(hSSC, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(hSSC$ssc_state), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839846 <- DimPlot(SRS9839846, cells = cells_to_use.SRS9839846, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839846)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839846$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839846@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839846.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839846<-NULL

# load query SRS9839847, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839847<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839847/R_SCT/SRS9839847.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839847) %in% rat_genes$Symbol
SRS9839847_v2 <- SRS9839847[tmp_idx,]
num<-length(rownames(SRS9839847_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839847_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839847_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839847_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839847_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839847_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839847<-SRS9839847_v2
DefaultAssay(SRS9839847) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839847, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839847 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839847, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839847$predicted.cluster <- factor(SRS9839847$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839847 <- DimPlot(SRS9839847, cells = cells_to_use.SRS9839847, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839847)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839847$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839847@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839847.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839847<-NULL

# load query SRS9839848, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839848<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839848/R_SCT/SRS9839848.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839848) %in% rat_genes$Symbol
SRS9839848_v2 <- SRS9839848[tmp_idx,]
num<-length(rownames(SRS9839848_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839848_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839848_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839848_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839848_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839848_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839848<-SRS9839848_v2
DefaultAssay(SRS9839848) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839848, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839848 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839848, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839848$predicted.cluster <- factor(SRS9839848$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839848 <- DimPlot(SRS9839848, cells = cells_to_use.SRS9839848, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839848)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839848$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839848@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839848.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839848<-NULL

# load query SRS9839849, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839849<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839849/R_SCT/SRS9839849.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839849) %in% rat_genes$Symbol
SRS9839849_v2 <- SRS9839849[tmp_idx,]
num<-length(rownames(SRS9839849_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839849_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839849_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839849_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839849_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839849_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839849<-SRS9839849_v2
DefaultAssay(SRS9839849) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839849, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839849 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839849, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839849$predicted.cluster <- factor(SRS9839849$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839849 <- DimPlot(SRS9839849, cells = cells_to_use.SRS9839849, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839849)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839849$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839849@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839849.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839849<-NULL

# load query SRS9839850, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839850<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839850/R_SCT/SRS9839850.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839850) %in% rat_genes$Symbol
SRS9839850_v2 <- SRS9839850[tmp_idx,]
num<-length(rownames(SRS9839850_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839850_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839850_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839850_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839850_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839850_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839850<-SRS9839850_v2
DefaultAssay(SRS9839850) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839850, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839850 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839850, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839850$predicted.cluster <- factor(SRS9839850$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839850 <- DimPlot(SRS9839850, cells = cells_to_use.SRS9839850, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839850)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839850$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839850@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839850.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839850<-NULL

# load query SRS9839851, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839851<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839851/R_SCT/SRS9839851.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839851) %in% rat_genes$Symbol
SRS9839851_v2 <- SRS9839851[tmp_idx,]
num<-length(rownames(SRS9839851_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839851_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839851_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839851_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839851_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839851_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839851<-SRS9839851_v2
DefaultAssay(SRS9839851) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839851, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839851 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839851, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839851$predicted.cluster <- factor(SRS9839851$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839851 <- DimPlot(SRS9839851, cells = cells_to_use.SRS9839851, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839851)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839851$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839851@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839851.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839851<-NULL

# load query SRS9839852, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839852<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839852/R_SCT/SRS9839852.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839852) %in% rat_genes$Symbol
SRS9839852_v2 <- SRS9839852[tmp_idx,]
num<-length(rownames(SRS9839852_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839852_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839852_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839852_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839852_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839852_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839852<-SRS9839852_v2
DefaultAssay(SRS9839852) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839852, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839852 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839852, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839852$predicted.cluster <- factor(SRS9839852$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839852 <- DimPlot(SRS9839852, cells = cells_to_use.SRS9839852, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839852)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839852$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839852@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839852.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839852<-NULL

# load query SRS9839853, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS9839853<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Rattus_norvegicus/SRS9839853/R_SCT/SRS9839853.rds')
rat_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/rat_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS9839853) %in% rat_genes$Symbol
SRS9839853_v2 <- SRS9839853[tmp_idx,]
num<-length(rownames(SRS9839853_v2))
genetable_1 <- data.frame(sSymbol=rownames(SRS9839853_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=rat_genes, by.x="sSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS9839853_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS9839853_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS9839853_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS9839853_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS9839853<-SRS9839853_v2
DefaultAssay(SRS9839853) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS9839853, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS9839853 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS9839853, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS9839853$predicted.cluster <- factor(SRS9839853$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS9839853 <- DimPlot(SRS9839853, cells = cells_to_use.SRS9839853, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("8-10 weeks\n(SRS9839853)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS9839853$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
write.table(SRS9839853@meta.data,file='states0to4.onto_which_projected_rat_germline.SRS9839853.metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
SRS9839853<-NULL

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/states0to4.onto_which_projected_rat_germline.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS9839846,fig.SRS9839847,fig.SRS9839848,fig.SRS9839849,fig.SRS9839850,fig.SRS9839851,fig.SRS9839852,fig.SRS9839853,ncol=3,nrow=3,common.legend=TRUE,legend='bottom')
dev.off()