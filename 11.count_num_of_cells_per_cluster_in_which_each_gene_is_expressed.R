## IMPORTANT: use "module add R-cbrg/202210" to ensure you use v4.1.1 of Seurat

### whole testes atlas ###

library(Seurat)

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 0)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster0.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster0.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster0.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 1)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 2)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 3)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 4)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster4.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster4.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster4.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 5)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster5.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster5.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster5.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 6)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster6.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster6.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster6.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 7)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster7.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster7.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster7.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 8)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster8.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster8.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster8.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

testes<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/human_adult.rds') # from 6.refine_clustering_of_whole_testes_atlas.R
cluster <- subset(testes, seurat_clusters == 9)
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='human_adult.number_of_cells_expressing_each_gene_in_Cluster9.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(cluster@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(cluster)
write.table(df,file='human_adult.proporp_of_cells_expressing_each_gene_in_Cluster9.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(cluster,assays='RNA')
write.table(df,file='human_adult.avg_expr_of_cells_expressing_each_gene_in_Cluster9.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

### SSC atlas ###

library(Seurat)

### resolution = 0.3 (4 clusters) ###

# cluster 0

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.0.3
ssc <- subset(states0to4, integrated_snn_res.0.3 == 0)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster0_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster0_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster0_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 1

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.0.3
ssc <- subset(states0to4, integrated_snn_res.0.3 == 1)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster1_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster1_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster1_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 2

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.0.3
ssc <- subset(states0to4, integrated_snn_res.0.3 == 2)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster2_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster2_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster2_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 3

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.0.3
ssc <- subset(states0to4, integrated_snn_res.0.3 == 3)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster3_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster3_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster3_resolution_0.3.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

### resolution = 1.1 (7 clusters) ###

# cluster 0

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 0)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster0_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster0_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster0_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 1

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 1)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster1_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster1_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster1_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 2

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 2)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster2_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster2_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster2_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 3

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 3)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster3_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster3_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster3_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 4

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 4)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster4_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster4_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster4_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 5

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 5)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster5_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster5_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster5_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 6

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
ssc <- subset(states0to4, integrated_snn_res.1.1 == 6)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster6_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster6_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster6_resolution_1.1.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

### resolution = 1.2 (8 clusters) ###

# cluster 0

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 0)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster0_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster0_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster0_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 1

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 1)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster1_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster1_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster1_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 2

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 2)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster2_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster2_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster2_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 3

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 3)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster3_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster3_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster3_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 4

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 4)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster4_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster4_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster4_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 5

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 5)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster5_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster5_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster5_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 6

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 6)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster6_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster6_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster6_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')

# cluster 7

states0to4<-readRDS('/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
Idents(states0to4)<-states0to4$integrated_snn_res.1.2
ssc <- subset(states0to4, integrated_snn_res.1.2 == 7)
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum')
write.table(df,file='states0to4.number_of_cells_expressing_each_gene_in_Cluster7_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-apply(ssc@assays$RNA@counts > 0, MARGIN = 1, FUN = 'sum') / ncol(ssc)
write.table(df,file='states0to4.proporp_of_cells_expressing_each_gene_in_Cluster7_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
df<-AverageExpression(ssc,assays='RNA')
write.table(df,file='states0to4.avg_expr_of_cells_expressing_each_gene_in_Cluster7_resolution_1.2.txt',row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')