## IMPORTANT: use "module add R-cbrg/202210" to ensure you use v4.1.1 of Seurat

library(Seurat)

states0to4<-readRDS('results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
DefaultAssay(states0to4) <- 'RNA' # it is important to set the DefaultAssay to RNA before running differential expression: "we don't recommend using the integrated matrix for differential expression" and "As a general rule, we always recommend performing DE on originally measured values - instead of on batch-corrected, imputed, etc. values. This ensures that the measurements that enter the DE test are indeed independent from each other, which is a requirement of any statistical DE test." (https://github.com/satijalab/seurat/issues/1057, https://github.com/satijalab/seurat/issues/1256 and https://github.com/satijalab/seurat/issues/2136)

##### ALL vs. ALL @ RESOLUTIONS 0.3, 1.1, AND 1.2 #####

### we wish to identify discriminative biomarkers of each cluster in the states 0-4 atlas, relative to all other clusters
### we first need to amend the identity class of states0to4 so that FindAllMarkers can be applied to it appropriately

Idents(states0to4)<-states0to4$integrated_snn_res.0.3
states0to4.markers <- FindAllMarkers(states0to4, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(states0to4.markers,file='states0to4.degs_All_vs_All_resolution_0.3.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
states0to4.markers <- FindAllMarkers(states0to4, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(states0to4.markers,file='states0to4.degs_All_vs_All_resolution_1.1.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
states0to4.markers <- FindAllMarkers(states0to4, test.use = "MAST", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(states0to4.markers,file='states0to4.degs_All_vs_All_resolution_1.2.txt',row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

##### X vs. Y and Y vs. X @ RESOLUTION 0.3 #####

### determine which genes are differentially expressed between the two subsets of undifferentiated SSCs in the r=0.3 atlas, in an X-vs-Y and Y-vs-X comparison
### these are clusters 1 and 2, respectively, 'A-dark' (EGR4+ PIWIL4+) and 'A-pale' (ID4+)

Idents(states0to4)<-states0to4$integrated_snn_res.0.3
deg <- FindMarkers(states0to4, ident.1 = "1", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_Adark_vs_Apale_resolution_0.3_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.0.3
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_Apale_vs_Adark_resolution_0.3_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

##### X vs. Y and Y vs. X @ RESOLUTION 1.1 #####

### determine which genes are differentially expressed between the three subsets of undifferentiated SSCs in the r=1.1 atlas, in an X-vs-Y and Y-vs-X comparison
### these are clusters 2, 0 and 4, respectively, 'state 0' (EGR4+ PIWIL4+), 'state 0A/1' (FGFR3+) and 'state 0B' (NANOS2+)

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_State0A_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "0", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_State0_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "4", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_State0B_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "4", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0B_vs_State0_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "0", ident.2 = "4", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_State0B_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "4", ident.2 = "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0B_vs_State0A_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

### determine which genes are differentially expressed between each of three subsets of undifferentiated SSCs *and* the subset of differentiating SSCs in the r=1.1 atlas, in an X-vs-Y and Y-vs-X comparison
### these are clusters 2, 0, 4 and 1, respectively, 'state 0' (EGR4+ PIWIL4+), 'state 0A/1' (FGFR3+) and 'state 0B' (NANOS2+) and 'SPG' (KIT+)
### rationale for this analysis: a differentiating cell is not a differentiated cell; it has not yet commited irreversibly to its fate and may still have the potential to do something else
### this explains an observation in mice: that cells sorted using the 'differentiating cell' marker MIWI2 (which is PIWIL4 in human) can be cultured in vitro, self-renewing as if they were 'undifferentiated stem cells' all along
### the goal of this DE analysis is to determine what genes are involved in a cell 'leaving the loop'

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "1", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_SPG_vs_State0_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_SPG_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "1", ident.2 = "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_SPG_vs_State0A_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "0", ident.2 = "1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_SPG_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "1", ident.2 = "4", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_SPG_vs_State0B_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.1
deg <- FindMarkers(states0to4, ident.1 = "4", ident.2 = "1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0B_vs_SPG_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

states0to4<-readRDS('results/states0to4.rds') # from 8.subset_whole_testes_atlas_to_obtain_States0to4.R
DefaultAssay(states0to4) <- 'RNA'
Idents(states0to4)<-states0to4$integrated_snn_res.1.1
new.cluster.ids <- c(
"state 0A/1",
"SPG",
"state 0",
"meiosis",
"state 0B",
"meiosis",
"meiosis"
)
names(new.cluster.ids) <- levels(states0to4)
states0to4 <- RenameIdents(states0to4, new.cluster.ids)

deg <- FindMarkers(states0to4, ident.1 = "state 0", ident.2 = "meiosis", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_meiosis_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "meiosis", ident.2 = "state 0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_meiosis_vs_State0_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "state 0A/1", ident.2 = "meiosis", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_meiosis_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "meiosis", ident.2 = "state 0A/1", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_meiosis_vs_State0A_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "state 0B", ident.2 = "meiosis", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0B_vs_meiosis_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "meiosis", ident.2 = "state 0B", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_meiosis_vs_State0B_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "SPG", ident.2 = "meiosis", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_SPG_vs_meiosis_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

deg <- FindMarkers(states0to4, ident.1 = "meiosis", ident.2 = "SPG", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_meiosis_vs_SPG_resolution_1.1_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

##### X vs. Y and Y vs. X @ RESOLUTION 1.2 #####

### determine which genes are differentially expressed between the four subsets of undifferentiated SSCs in the r=1.2 atlas, in an X-vs-Y and Y-vs-X comparison
### these are clusters 2, 5, 0 and 6
### 2 = 'state 0' (EGR4+ PIWIL4+), 5 = 'state 0/0B', 6 = 'state 0B/1', 0 = 'state 0A/1' (FGFR3+)

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_State0A_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "5", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_State0to0B_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "2", ident.2 = "6", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0_vs_State0Bto1_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "0", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_State0_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "0", ident.2 = "5", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_State0to0B_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "0", ident.2 = "6", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0A_vs_State0Bto1_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "5", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0to0B_vs_State0_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "5", ident.2 = "6", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0to0B_vs_State0Bto1_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "5", ident.2 = "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0to0B_vs_State0A_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "6", ident.2 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0Bto1_vs_State0_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "6", ident.2 = "5", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0Bto1_vs_State0to0B_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')

Idents(states0to4)<-states0to4$integrated_snn_res.1.2
deg <- FindMarkers(states0to4, ident.1 = "6", ident.2 = "0", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(deg,file='states0to4.degs_State0Bto1_vs_State0A_resolution_1.2_onlyPosTRUE_minPCT25_logFC25.txt',row.names=TRUE,col.names=TRUE,quote=FALSE,sep='\t')