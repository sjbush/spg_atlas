## PROJECTIONS: of three FACS-sorted steady state spermatogonial samples onto the states 0-4 UMAP. This is a validation test of the gating strategy to distinguish somatic from germline cells, before re-projecting the latter onto the states 0-4 UMAP
## IMPORTANT: use "module add R-cbrg/202210" to ensure you use v4.1.1 of Seurat
## WARNING: "there are multiple ggarranges. I was trying to use the command with the egg package, but common.legend only works with the one in ggpubr" https://stackoverflow.com/questions/67305903/no-applicable-method-for-ggplot-build-applied-to-an-object-of-class-logical

library(Seurat)
library(ggpubr)
library(tidyverse)
theme_set(theme_bw())

# STEP 1: project FACS-sorted Hermann spermatagonial samples onto the adult whole testes UMAP

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

# load query SRS2823404 and then project onto the adult whole testes UMAP
SRS2823404<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823404/R_SCT/SRS2823404.rds')
DefaultAssay(SRS2823404) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS2823404,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS2823404, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS2823404 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS2823404, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("SSC","spermatogonia","spermatocyte","early spermatid (1)","early spermatid (2)","late spermatid (1)","late spermatid (2)","myoid/Leydig","Sertoli","endothelia") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS2823404$predicted.cluster <- factor(SRS2823404$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(adult, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(adult$ssc_state))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS2823404 <- DimPlot(SRS2823404, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("SRS2823404") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS2823404$predicted.cluster))

# load query SRS2823405 and then project onto the adult whole testes UMAP
SRS2823405<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823405/R_SCT/SRS2823405.rds')
DefaultAssay(SRS2823405) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS2823405,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS2823405, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS2823405 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS2823405, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("SSC","spermatogonia","spermatocyte","early spermatid (1)","early spermatid (2)","late spermatid (1)","late spermatid (2)","myoid/Leydig","Sertoli","endothelia") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS2823405$predicted.cluster <- factor(SRS2823405$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS2823405 <- DimPlot(SRS2823405, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("SRS2823405") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS2823405$predicted.cluster))

# load query SRS2823406 and then project onto the adult whole testes UMAP
SRS2823406<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823406/R_SCT/SRS2823406.rds')
DefaultAssay(SRS2823406) <- "SCT"
DefaultAssay(adult) <- "integrated"
adult@meta.data$ssc_state <- as.character(adult@meta.data$ssc_state)
common_features <- lapply(list(SRS2823406,adult), row.names) %>% Reduce(intersect, .)
adult.anchors <- FindTransferAnchors(reference = adult, query = SRS2823406, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
adult <- RunUMAP(adult, dims = 1:prin_comp, reduction = "pca", n.neighbors = 20, return.model = TRUE)
SRS2823406 <- MapQuery(anchorset = adult.anchors, reference = adult, query = SRS2823406, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("SSC","spermatogonia","spermatocyte","early spermatid (1)","early spermatid (2)","late spermatid (1)","late spermatid (2)","myoid/Leydig","Sertoli","endothelia") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
adult$ssc_state <- factor(adult$ssc_state, levels = levels.set)
SRS2823406$predicted.cluster <- factor(SRS2823406$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS2823406 <- DimPlot(SRS2823406, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("SRS2823406") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-15, 15)) + scale_y_continuous(name="UMAP_2", limits=c(-15, 15)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_discrete(drop=TRUE,limits = levels(SRS2823406$predicted.cluster))

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_Hermann_samples.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS2823404,fig.SRS2823405,fig.SRS2823406,ncol=2,nrow=2,common.legend=TRUE,legend='right')
dev.off()

fig1<-ggarrange(fig.ref,fig.SRS2823404,fig.SRS2823405,fig.SRS2823406,ncol=2,nrow=2,common.legend=TRUE,legend='right')
#fig1<-annotate_figure(plot, top = text_grob("A", color = "red", face = "bold", size = 14))

predicted_sscs.SRS2823404 <- subset(SRS2823404@meta.data,((SRS2823404@meta.data$predicted.cluster != "myoid/Leydig") & (SRS2823404@meta.data$predicted.cluster != "Sertoli") & (SRS2823404@meta.data$predicted.cluster != "endothelia")))
predicted_sscs.SRS2823405 <- subset(SRS2823405@meta.data,((SRS2823405@meta.data$predicted.cluster != "myoid/Leydig") & (SRS2823405@meta.data$predicted.cluster != "Sertoli") * (SRS2823405@meta.data$predicted.cluster != "endothelia")))
predicted_sscs.SRS2823406 <- subset(SRS2823406@meta.data,((SRS2823406@meta.data$predicted.cluster != "myoid/Leydig") & (SRS2823406@meta.data$predicted.cluster != "Sertoli") & (SRS2823406@meta.data$predicted.cluster != "endothelia")))
write.table(predicted_sscs.SRS2823404,file='/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_SRS2823404.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
write.table(predicted_sscs.SRS2823405,file='/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_SRS2823405.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
write.table(predicted_sscs.SRS2823406,file='/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_SRS2823406.cells_which_project_onto_germline.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

# STEP 2: project FACS-sorted Hermann spermatagonial samples onto the states 0-4 UMAP, after first having projected the same sample to the whole-adult atlas so that we can exclude those cells which do NOT project to the SSC compartment

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
# we had subset that list to contain only those cells which map to the germline clusters of the whole-adult atlas
# we will now constrain all subsequent (re)projections onto to these cells
# the way to do this is NOT to subset the Seurat object as that will leave you with too few cells (the k.score and k.filter parameters of FindTransferAnchors will cause errors); rather, we tell DimPlot to only print certain cells
predicted_sscs.SRS2823404<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_SRS2823404.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS2823405<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_SRS2823405.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')
predicted_sscs.SRS2823406<-read.table('/project/GorielyLab2021/sbush/human_adult_SSCs/results/adult_whole_testes.onto_which_projected_SRS2823406.cells_which_project_onto_germline.txt',header=TRUE,sep='\t')

cells_to_use.SRS2823404<-predicted_sscs.SRS2823404$barcode
cells_to_use.SRS2823405<-predicted_sscs.SRS2823405$barcode
cells_to_use.SRS2823406<-predicted_sscs.SRS2823406$barcode

# load query SRS2823404 and then project onto the SSC state 0-4 UMAP
SRS2823404<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823404/R_SCT/SRS2823404.rds')
DefaultAssay(SRS2823404) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
common_features <- lapply(list(SRS2823404,hSSC), row.names) %>% Reduce(intersect, .)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS2823404, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS2823404 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS2823404, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS2823404$predicted.cluster <- factor(SRS2823404$predicted.cluster, levels = levels.set)
fig.ref <- DimPlot(hSSC, reduction = "umap", group.by = "ssc_state", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("Reference") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(hSSC$ssc_state), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS2823404 <- DimPlot(SRS2823404, cells = cells_to_use.SRS2823404, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("SRS2823404") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS2823404$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))

# load query SRS2823405 and then project onto the SSC state 0-4 UMAP
SRS2823405<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823405/R_SCT/SRS2823405.rds')
DefaultAssay(SRS2823405) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
common_features <- lapply(list(SRS2823405,hSSC), row.names) %>% Reduce(intersect, .)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS2823405, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS2823405 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS2823405, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS2823405$predicted.cluster <- factor(SRS2823405$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS2823405 <- DimPlot(SRS2823405, cells = cells_to_use.SRS2823405, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("SRS2823405") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS2823405$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))

# load query SRS2823406 and then project onto the SSC state 0-4 UMAP
SRS2823406<-readRDS('/project/GorielyLab2021/sbush/ssc_atlas/kb/Homo_sapiens/SRS2823406/R_SCT/SRS2823406.rds')
DefaultAssay(SRS2823406) <- "SCT"
DefaultAssay(hSSC) <- "integrated"
hSSC@meta.data$ssc_state <- as.character(hSSC@meta.data$ssc_state)
common_features <- lapply(list(SRS2823406,hSSC), row.names) %>% Reduce(intersect, .)
hSSC.anchors <- FindTransferAnchors(reference = hSSC, query = SRS2823406, dims = 1:prin_comp, reference.reduction = "pca", normalization.method="SCT", features=common_features)
hSSC <- RunUMAP(hSSC, dims = 1:prin_comp, reduction = "pca", n.neighbors = 200, return.model = TRUE)
SRS2823406 <- MapQuery(anchorset = hSSC.anchors, reference = hSSC, query = SRS2823406, refdata = list(cluster = "ssc_state"), reference.reduction = "pca", reduction.model = "umap")
levels.set = c("state 0","state 0A/1","state 0B","SPG","diff SPG","leptotene","zygotene") # we need to set "ssc_state" and "predicted.cluster" as two factors with the same level. This way, the two UMAPs, p1 and p2, will show matched colours. See https://github.com/satijalab/seurat/issues/5582
hSSC$ssc_state <- factor(hSSC$ssc_state, levels = levels.set)
SRS2823406$predicted.cluster <- factor(SRS2823406$predicted.cluster, levels = levels.set)
# we use 'limits' to ensure we plot exactly the same colour for each level of the factors in $ssc_state and $predicted.cluster, even when not every level appears in both UMAPs. See https://stackoverflow.com/questions/6919025/how-to-assign-colors-to-categorical-variables-in-ggplot2-that-have-stable-mappin
fig.SRS2823406 <- DimPlot(SRS2823406, cells = cells_to_use.SRS2823406, reduction = "ref.umap", group.by = "predicted.cluster", pt.size = 0.2, label = FALSE) + NoLegend() + ggtitle("SRS2823406") + theme(aspect.ratio=1) + scale_x_continuous(name="UMAP_1", limits=c(-12, 8)) + scale_y_continuous(name="UMAP_2", limits=c(-10, 10)) + theme(plot.title = element_text(size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 8)) + scale_colour_manual(drop=TRUE,limits = levels(SRS2823406$predicted.cluster), values=c("state 0" = "orange", "state 0A/1" = "green", "state 0B" = "red", "SPG" = "yellow", "diff SPG" = "pink", "leptotene" = "cornflowerblue", "zygotene" = "blueviolet"))

pdf('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.onto_which_projected_Hermann_samples_after_filtering.pdf',onefile=FALSE)
ggarrange(fig.ref,fig.SRS2823404,fig.SRS2823405,fig.SRS2823406,ncol=2,nrow=2,common.legend=TRUE,legend='right')
dev.off()

fig2<-ggarrange(fig.ref,fig.SRS2823404,fig.SRS2823405,fig.SRS2823406,ncol=2,nrow=2,common.legend=TRUE,legend='right')
#fig2<-annotate_figure(plot, top = text_grob("B", color = "red", face = "bold", size = 14))

# STEP 3: print a combined figure

fig<-ggarrange(fig1, fig2, nrow=1, ncol=2)
ggsave('figure_s12.pdf',plot=fig,width=14,height=7)