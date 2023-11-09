# code from https://www.kallistobus.tools/tutorials/kb_building_atlas/r/kb_analysis_0_r
# for published write-ups of this methodology see https://www.cell.com/cell-reports-medicine/fulltext/S2666-3791(21)00253-6 (Di Persio 2021), https://www.nature.com/articles/s41422-018-0099-2#Sec19 (Guo 2018) and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6402825/ (Sohni 2019)
# see also the Seurat cheat sheet: https://satijalab.org/seurat/articles/essential_commands.html
# for an overview, see: https://hbctraining.github.io/scRNA-seq/lessons/02_SC_generation_of_count_matrix.html

# IMPORTANT: HARD-CODED ELEMENTS OF THIS SCRIPT
# the t2g file is at /project/GorielyLab2021/sbush/ssc_atlas/indexes/Ovis_aries/t2g2.txt

library(Seurat)
library(Matrix)
library(glmGamPoi)
library(tidyverse)
library(patchwork)
library(sctransform)
library(clustree)
theme_set(theme_bw())

read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

get_knee_df <- function(mat) {
  total <- rank <- NULL
  tibble(total = Matrix::colSums(mat),
         rank = row_number(desc(total))) %>%
    distinct() %>%
    dplyr::filter(total > 0) %>% 
    arrange(rank)
}

get_inflection <- function(df, lower = 100) {
  log_total <- log_rank <- total <-  NULL
  df_fit <- df %>% 
    dplyr::filter(total > lower) %>% 
    transmute(log_total = log10(total),
              log_rank = log10(rank))
  d1n <- diff(df_fit$log_total)/diff(df_fit$log_rank)
  right.edge <- which.min(d1n)
  10^(df_fit$log_total[right.edge])
}

knee_plot <- function(df, inflection) {
total <- rank_cutoff <- NULL
annot <- tibble(inflection = inflection, rank_cutoff = max(df$rank[df$total > inflection]))
ggplot(df, aes(total, rank)) +
    geom_path() +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2, color = "gray40") +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2,color = "gray40") +
    geom_text(aes(inflection, rank_cutoff, label = paste(rank_cutoff, "'cells'")), data = annot, vjust = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks()
}

### read the directory of UNFILTERED data

res_mat <- read_count_output("./counts_unfiltered", name = "cells_x_genes")
num_cells_in_kb_unfiltered <- ncol(res_mat)

### OUT1: make a knee plot from the UNFILTERED data

pdf('R_SCT/knee_plot_unfiltered.pdf')
options(repr.plot.width=9, repr.plot.height=6)
knee_df <- get_knee_df(res_mat)
inflection <- get_inflection(knee_df)
knee_plot(knee_df, inflection)
dev.off()

### OUT2: output a table stating how many cells have x UMIs: this is the basis by which KB implements its filter

num_umis_at_inflection_point<-as.numeric(inflection[[1]])
num_cells_at_inflection_point<-max(knee_df$rank[knee_df$total > inflection])
inflection_df<-data.frame(UMIS=num_umis_at_inflection_point,CELLS=num_cells_at_inflection_point)
write.table(inflection_df,file='R_SCT/inflection_point.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

### read the directory of FILTERED data

res_mat <- read_count_output("./counts_filtered", name = "cells_x_genes")
num_cells_in_kb_filtered <- ncol(res_mat)
num_genes_in_kb_filtered <- nrow(res_mat)

### FILTER CELLS TO REMOVE UNWANTED ONES FROM THE DATASET. WE WILL EMPIRICALLY DERIVE APPROPRIATE VALUES WITH THE FOLLOWING TWO FIGURES.
# for definitions, see https://www.biostars.org/p/407036/. nFeature_RNA = number of genes detected in each cell. nCount_RNA = total number of molecules detected within a cell.
# how many cells do we have left? See testis@meta.data, which we will export below
# N.B. example parameter choices of nFeature_RNA (number of detected features [genes] per cell: 500 to 10000), nCount_RNA (number of UMIs expressed per cell: 1000 to 50000), and percentage of mitochondrial genes expressed (<0.2%) are derived from https://www.pnas.org/content/pnas/suppl/2020/07/13/2000362117.DCSupplemental/pnas.2000362117.sapp.pdf

tr2g <- read_tsv("/project/GorielyLab2021/sbush/ssc_atlas/indexes/Ovis_aries/t2g2.txt", col_names = c("transcript", "gene", "gene_name")) # HARD-CODED
tr2g <- distinct(tr2g[, c("gene", "gene_name")])
rownames(res_mat) <- tr2g$gene_name[match(rownames(res_mat), tr2g$gene)]
(testis <- CreateSeuratObject(counts = res_mat, project = "testis", min.cells = 3, min.features = 100))
testis[["percent.mt"]] <- PercentageFeatureSet(testis, pattern = "^MT-")

### OUT3: visualise QC metrics as a violin plot

pdf('R_SCT/violin_plots_of_QC_metrics.pdf')
options(repr.plot.width=12, repr.plot.height=6)
VlnPlot(testis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

### OUT4: create scatter plots of the % of MT genes, no. of UMIs in each cell, and no. of distinct genes in each cell (multiple UMIs - i.e. transcripts - are available per gene)

# we will now filter cells on the basis of the proportion of mitochondrial genes and the number of RNA molecules captured, excluding those with UMI counts that are too low (indicative of empty or low-quality droplets) or too high (indicative of doublets)
# to do this, we'll generate some figures in order to work out appropriate values

pdf('R_SCT/scatter_plots_to_inform_filtering.pdf')
plot1 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
plot1 + plot2
dev.off()

### APPLY SCTRANSFORM NORMALIZATION, REGRESSING OUT THE EFFECT OF SEVERAL VARIABLES

testis <- SCTransform(testis, method = "glmGamPoi", vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt")) # from https://github.com/satijalab/seurat/issues/1739: it is not really necessary to set UMI counts ("nCount_RNA") in the vars.to.regress because these are already set in SCTransform's regularized negative binomial model

# we will also mitigate the effect of cell-cycle heterogeneity by calculating cell cycle phase scores so that we can regress these out too
# this process is described in https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html and https://github.com/satijalab/seurat/issues/1679, with https://www.biostars.org/p/364214/ detailing the location of lists of cell cycle genes. Some are from Table S2 of Dominguez 2018 (https://www.nature.com/articles/cr201684) but others are pre-loaded in Seurat
# note potential issues with applying the CellCycleScoring function in conjunction with SCTransform (namely, with small datasets, the error message "insufficient data values to produce 24 bins"): https://github.com/satijalab/seurat/issues/1227 and https://github.com/satijalab/seurat/issues/3692
s.genes <- c("E2F1","CCNE1","CCNE2","POLD3","DTL","MCM6","PCNA","CHAF1A","SLBP","UNG","CDC6","MCM2","PASK","HSPB8","EXO1","RFC4") # from Dominguez 2018. N.B. Seurat provides its own list at cc.genes.updated.2019$s.genes
g2m.genes <- c("GTSE1","GPSM2","CCNF","KIF11","H2AFX","CDCA8","NDE1","ESPL1","CCNA2","KIF23","GPR126","KIF22","ARHGAP11A","BUB3","HMGB2","TOP2A","BRD8","UBE2C","NUSAP1","CKS1B","CDC25C","FAM64A","CENPE","CKAP2","NEK2","CDCA3","MKI67","BUB1","TPX2","PLK1","CKS2","UBE2S","CDC20","CCNB1","KPNA2","TTK","HMMR","SPAG5","CENPF","CDC25B","HMGB3","KIF2C","TACC3","BUB1B","CCNB2","LBR","BIRC5","PTTG1","TROAP","PRC1","SFPQ") # from Dominguez 2018. N.B. Seurat provides its own list at g2m.genes <- cc.genes.updated.2019$g2m.genes
testis <- NormalizeData(testis, assay = 'RNA')
testis <- CellCycleScoring(testis, s.features = s.genes, g2m.features = g2m.genes, assay = 'RNA', set.ident = TRUE)
testis <- SCTransform(testis, method = "glmGamPoi", assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score"))

### OUT5: determine the dimensionality of the data (1). Perform PCA analysis on the scaled data and generate an elbow plot: a ranking of principal components by the percentage of variance explained by each ones

testis <- RunPCA(testis, verbose = FALSE) # note that the default number of principal components analysed is 50 and that there must be more principal components than there are cells. See https://github.com/satijalab/seurat/issues/1914

pdf('R_SCT/elbow_plot.pdf')
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(testis, ndims = 50)
dev.off()

# quantify content of the elbow plot programmatically. To do this, we implement code from https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html

pct <- testis[["pca"]]@stdev / sum(testis[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%

# let us choose the minimum of these two metrics and conclude that at this point the PCs cover the majority of the variation in the data. We will use all PCs up to and including this 'elbow' for downstream analyses, i.e. dimension reduction and UMAP visualisation.
prin_comp <- min(component1, component2)
write.table(prin_comp,file='R_SCT/elbow_PC.txt',row.names=FALSE,col.names=FALSE,quote=FALSE,sep='\t')

### OUT6: VISUALISE THE CLUSTERING WITH UMAP, USING AS INPUT THE DIMENSIONALITY INFERRED FROM THE ELBOW PLOT

# see description of UMAP hyperparameters here: https://smorabit.github.io/blog/2020/umap/
# when it comes to clustering the UMAP, we're going to vary the resolution parameter of the FindClusters function in order to later optimise its value using clustree, as explained here: https://github.com/lazappi/clustree/issues/18
# for additional detail on clustree, see also http://oshlacklab.com/combes-organoid-paper/04_Organoids_Clustering.html, https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#clustering-trees-for-scrna-seq-data, and the clustree paper, Zappia and Oshlack (2018): https://academic.oup.com/gigascience/article/7/7/giy083/5052205?login=true

resolution.range <- seq(from = 0, to = 2.0, by = 0.2)

pdf('R_SCT/umap_plot.pdf')
testis <- FindNeighbors(testis, reduction = 'pca', dims = 1:prin_comp, verbose = FALSE)
testis <- FindClusters(testis, algorithm=3, resolution = resolution.range, verbose = FALSE) # from https://satijalab.org/seurat/archive/v3.1/pbmc3k_tutorial.html: the FindClusters function "contains a resolution parameter that sets the '˜granularity' of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets"
testis <- RunUMAP(testis, dims = 1:prin_comp, verbose = FALSE) # note that the default values for k.param in FindNeighbors (20) and n.neighbors in RunUMAP (30) are different, but both define the number of neighbouring points. Discussed at https://github.com/satijalab/seurat/issues/4717: "if you want to ensure that the UMAP representation and clustering results are as consistent as possible it's a good idea to use the same number of nearest neighbors when building the neighbor graph for clustering and for UMAP". k.param is an especially important parameter with significant effect upon clustering: https://jtggjournal.com/article/view/3809
tplot = DimPlot(testis, reduction = "umap", label=TRUE, pt.size = .1)
tplot[[1]]$layers[[1]]$aes_params$alpha = 0.5
print(tplot)
dev.off()

### OUT7: export the Clustree plot

# how to interpret a Clustree plot? See figure 4 of https://academic.oup.com/gigascience/article/7/7/giy083/5052205 and the links above
# the size of each node (i.e. cluster) is related to the number of cells it contains.
# the colour of each node indicates the clustering resolution
# nodes are named according to the Seurat cluster ID, assigned by the UMAP, i.e. the contents of SRS6959440@meta.data$seurat_clusters. If you read the plot horizontally, you can see that for each of the different resolutions (node colours), there are a different set of cluster names. At resolution = 0, there is one cluster, number 0. At resolution = 0.1, there are 8 clusters, numbered 0 through 7, and so on.
# edges are coloured according to the incoming node proportion, i.e. the number of cells in the edge divided by the number of cells in the node it points to. This metric shows the importance of the edge to the higher-resolution cluster independently of the cluster size.
# at increasingly high resolutions, we expect to see the tree becoming messier, containing nodes with multiple incoming edges *and* an increasing number of nodes with *low* in-proportion edges, together indicating cluster instability. What this means is that as we increase resolution, smaller-sized clusters will become increasingly fragmented as a modest proportion of their cells will be partitioned into separate clusters.
# we can also interpret the Clustree plot using the 'SC3 stability score': "starting with a set of cluster labels at different resolutions, each cluster is scored, with clusters awarded increased stability if they share the same samples as a cluster at another resolution but penalized for being at a higher resolution" (see https://academic.oup.com/gigascience/article/7/7/giy083/5052205)

pdf('R_SCT/clustree_plot.pdf',paper="a4r")
clustree_fig<-clustree(testis,prefix="SCT_snn_res.")
print(clustree_fig)
dev.off()

### OUT8: export the Seurat metadata

# note that this table contains the barcode for each cell. We can subset the Seurat object to exclude cells by barcode: https://github.com/satijalab/seurat/issues/1693
# the number of rows in this table = the number of cells we have used

barcodes<-rownames(testis@meta.data)
testis@meta.data$barcode<-barcodes
write.table(testis@meta.data,file='R_SCT/seurat_metadata.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

### OUT9: output a table showing how many cells and how many genes we have

total_genes<-nrow(x=testis)
total_cells<-ncol(x=testis) # = num_cells_passing_all_QC_thresholds
total_df<-data.frame(NUM_CELLS_KB_UNFILTERED=num_cells_in_kb_unfiltered,NUM_CELLS_KB_FILTERED=num_cells_in_kb_filtered,NUM_GENES_KB_FILTERED=num_genes_in_kb_filtered,NUM_CELLS_WITH_GE_100_GENES=total_cells,NUM_GENES_IN_GE_3_CELLS=total_genes)
write.table(total_df,file='R_SCT/totals.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

### OUT10: export the Seurat object for later integration

saveRDS(testis, file = "R_SCT/testis.rds")