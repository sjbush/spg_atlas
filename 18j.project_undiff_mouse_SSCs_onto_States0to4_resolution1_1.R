## PROJECTIONS
## IMPORTANT: use "module add R-cbrg/202210" to ensure you use v4.1.1 of Seurat
## NOTE: Hobbs data is from La 2018 "Identification of dynamic undifferentiated cell states within the male germline" https://www.nature.com/articles/s41467-018-04827-z (GEO: GSE107256 | ENA: PRJNA449375)
## Description of this dataset from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA449375: "Undifferentiated cells were isolated from adult mouse testes based on expression of a Plzf/Zbtb16 gene reporter, which marks both stem and progenitor spermatogonia"

# STEP 1: project the two Hobbs undifferentiated mouse SSC samples onto the adult whole testes UMAP

library(Seurat)
library(ggpubr)
library(tidyverse)
theme_set(theme_bw())

# load reference: adult human whole testes data
adult <- readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
new.cluster.ids <- c(
"late spermatid (1)",
"early spermatid (1)",
"early spermatid (2)",
"myoid/Leydig",
"spermatocyte",
"late spermatid (2)",
"SSC",
"endothelia",
"spermatogonia",
"Sertoli"
)
new.cluster.ids <- factor(new.cluster.ids, levels = c("SSC","spermatogonia","spermatocyte","early spermatid (1)","early spermatid (2)","late spermatid (1)","late spermatid (2)","myoid/Leydig","Sertoli","endothelia"))
names(new.cluster.ids) <- levels(adult)
adult <- RenameIdents(adult, new.cluster.ids)
adult@meta.data$ssc_state <- Idents(adult)

# to re-create the human UMAP, we first need to quantify the content of the elbow plot, implementing code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- adult[["pca"]]@stdev / sum(adult[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let us choose the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data. We need this parameter to input into RunUMAP
prin_comp <- min(component1, component2)

# load query SRS3150566, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3150566<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3150566/R_SCT/SRS3150566.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3150566) %in% mouse_genes$Symbol
SRS3150566_v2 <- SRS3150566[tmp_idx,]
num<-length(rownames(SRS3150566_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3150566_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3150566_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3150566_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3150566_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3150566_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3150566<-SRS3150566_v2
DefaultAssay(SRS3150566) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3150566, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3150566 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3150566, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("SSC","spermatogonia","spermatocyte","early spermatid (1)","early spermatid (2)","late spermatid (1)","late spermatid (2)","myoid/Leydig","Sertoli","endothelia") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3150566$predicted.cluster <- factor(SRS3150566$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.1, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 10)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3150566 <- DimPlot(SRS3150566, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.1, label = FALSE) + NoLegend() + ggtitle("undifferentiated mouse SSCs\n(SRS3150566)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 10)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3150566$predicted.cluster))

# load query SRS3150567, edit gene names to match that of human, and then project onto the SSC state 0-4 UMAP
SRS3150567<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3150567/R_SCT/SRS3150567.rds')
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3150567) %in% mouse_genes$Symbol
SRS3150567_v2 <- SRS3150567[tmp_idx,]
num<-length(rownames(SRS3150567_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3150567_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3150567_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3150567_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3150567_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3150567_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3150567<-SRS3150567_v2
DefaultAssay(SRS3150567) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS3150567, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT")
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS3150567 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS3150567, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("SSC","spermatogonia","spermatocyte","early spermatid (1)","early spermatid (2)","late spermatid (1)","late spermatid (2)","myoid/Leydig","Sertoli","endothelia") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS3150567$predicted.cluster <- factor(SRS3150567$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3150567 <- DimPlot(SRS3150567, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.1, label = FALSE) + NoLegend() + ggtitle("undifferentiated mouse SSCs\n(SRS3150567)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 10)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS3150567$predicted.cluster))

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_undiff_mouse_SSCs.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS3150566,fig.SRS3150567,nrow=2,ncol=2,common.legend=TRUE,legend='right')
dev.off()

fig1<-ggarrange(fig.ref,fig.SRS3150566,fig.SRS3150567,nrow=2,ncol=2,common.legend=TRUE,legend='right')

predicted_sscs.SRS3150566 <- subset(SRS3150566@meta.data,((SRS3150566@meta.data$predicted.cluster != "myoid/Leydig") & (SRS3150566@meta.data$predicted.cluster != "Sertoli") & (SRS3150566@meta.data$predicted.cluster != "endothelia")))
predicted_sscs.SRS3150567 <- subset(SRS3150567@meta.data,((SRS3150567@meta.data$predicted.cluster != "myoid/Leydig") & (SRS3150567@meta.data$predicted.cluster != "Sertoli") & (SRS3150567@meta.data$predicted.cluster != "endothelia")))
write.table(predicted_sscs.SRS3150566,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_SRS3150566.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
write.table(predicted_sscs.SRS3150567,file='/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_SRS3150567.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

######################################################################

# STEP 2: project the two Hobbs undifferentiated mouse SSC samples onto the adult states 0-4 UMAP, after first having projected the same sample to the whole-adult atlas so that we can exclude those cells which do NOT project to the SSC compartment

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

# we have previously obtained the list of predicted cluster identities after projecting each query pre-pubertal sample onto the whole-adult atlas (which has been clustered at resolution 0.15)
# we had subset that list to contain only those cells which map to whole-adult atlas cluster 6, the SSCs
# we will now constrain all subsequent (re)projections to these cells
# the way to do this is NOT to subset the Seurat object as that will leave you with too few cells (the k.score and k.filter parameters of FindTransferAnchors will cause errors); rather, we tell DimPlot to only print certain cells
predicted_sscs.SRS3150566<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_SRS3150566.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS3150567<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/adult_whole_testes.onto_which_projected_SRS3150567.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')

cells_to_use.SRS3150566<-predicted_sscs.SRS3150566$barcode
cells_to_use.SRS3150567<-predicted_sscs.SRS3150567$barcode

# load query SRS3150566 and then project onto the SSC state 0-4 UMAP
SRS3150566 <- readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3150566/R_SCT/SRS3150566.rds')
SRS3150566 <- FindNeighbors(SRS3150566, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
SRS3150566 <- FindClusters(SRS3150566, algorithm=3, resolution = 0.2, verbose = FALSE)
SRS3150566 <- RunUMAP(SRS3150566, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3150566) %in% mouse_genes$Symbol
SRS3150566_v2 <- SRS3150566[tmp_idx,]
num<-length(rownames(SRS3150566_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3150566_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3150566_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3150566_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3150566_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3150566_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3150566<-SRS3150566_v2
DefaultAssay(SRS3150566) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
common_features <- lapply(list(SRS3150566,hSSC), row.names) %>% Reduce(intersect, .)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS3150566, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS3150566 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS3150566, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS3150566$predicted.cluster <- factor(SRS3150566$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(hSSC, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(hSSC$ssc_state), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3150566 <- DimPlot(SRS3150566, cells = cells_to_use.SRS3150566, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("undifferentiated mouse SSCs\n(SRS3150566)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS3150566$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))

# load query SRS3150567 and then project onto the SSC state 0-4 UMAP
SRS3150567 <- readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Mus_musculus/SRS3150567/R_SCT/SRS3150567.rds')
SRS3150567 <- FindNeighbors(SRS3150567, reduction = 'pca', dims = 1:prin_comp, k.param = 20, verbose = FALSE)
SRS3150567 <- FindClusters(SRS3150567, algorithm=3, resolution = 0.2, verbose = FALSE)
SRS3150567 <- RunUMAP(SRS3150567, dims = 1:prin_comp, n.neighbors = 20, verbose = FALSE)
mouse_genes<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/gene_lists/mouse_to_human_symbol_lookup.txt',header=T,sep='\t')
tmp_idx <- rownames(SRS3150567) %in% mouse_genes$Symbol
SRS3150567_v2 <- SRS3150567[tmp_idx,]
num<-length(rownames(SRS3150567_v2))
genetable_1 <- data.frame(mSymbol=rownames(SRS3150567_v2), order=c(1:num), stringsAsFactors = F)
genetable_2 <- merge(x=genetable_1, y=mouse_genes, by.x="mSymbol", by.y="Symbol")
genetable_2 <- genetable_2[order( genetable_2$order ), ]
rownames(SRS3150567_v2@assays$RNA@counts) <- genetable_2$HumanSymbol
rownames(SRS3150567_v2@assays$RNA@data) <- genetable_2$HumanSymbol
rownames(SRS3150567_v2@assays$SCT@counts) <- genetable_2$HumanSymbol
rownames(SRS3150567_v2@assays$SCT@data) <- genetable_2$HumanSymbol
SRS3150567<-SRS3150567_v2
DefaultAssay(SRS3150567) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
common_features <- lapply(list(SRS3150567,hSSC), row.names) %>% Reduce(intersect, .)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS3150567, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS3150567 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS3150567, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS3150567$predicted.cluster <- factor(SRS3150567$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS3150567 <- DimPlot(SRS3150567, cells = cells_to_use.SRS3150567, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("undifferentiated mouse SSCs\n(SRS3150567)") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS3150567$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/states0to4.onto_which_projected_undiff_mouse_SSCs.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS3150566,fig.SRS3150567,nrow=2,ncol=2,common.legend=TRUE,legend='bottom')
dev.off()

fig2<-ggarrange(fig.ref,fig.SRS3150566,fig.SRS3150567,nrow=2,ncol=2,common.legend=TRUE,legend='bottom')

######################################################################

# STEP 3: print a combined figure

fig<-ggarrange(fig1, fig2, nrow=1, ncol=2)
ggsave('figure_s13.pdf',plot=fig,width=14,height=7)