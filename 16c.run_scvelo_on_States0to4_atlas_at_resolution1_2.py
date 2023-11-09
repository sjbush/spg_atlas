import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

sample1_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS6959440.csv")
sample1 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS6959440/counts_unfiltered/adata.loom")
sample1 = sample1[np.isin(sample1.obs.barcode,sample1_barcodes["x"])]
sample1.obs['barcodeID'] = sample1.obs.barcode
for i in range(len(sample1.obs.barcode)):
   sample1.obs['barcodeID'][i] = sample1.obs.barcode[i] + '_1'

sample2_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS6959441.csv")
sample2 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS6959441/counts_unfiltered/adata.loom")
sample2 = sample2[np.isin(sample2.obs.barcode,sample2_barcodes["x"])]
sample2.obs['barcodeID'] = sample2.obs.barcode
for i in range(len(sample2.obs.barcode)):
   sample2.obs['barcodeID'][i] = sample2.obs.barcode[i] + '_2'

sample3_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS6959442.csv")
sample3 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS6959442/counts_unfiltered/adata.loom")
sample3 = sample3[np.isin(sample3.obs.barcode,sample3_barcodes["x"])]
sample3.obs['barcodeID'] = sample3.obs.barcode
for i in range(len(sample3.obs.barcode)):
   sample3.obs['barcodeID'][i] = sample3.obs.barcode[i] + '_3'

sample4_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS4181127.csv")
sample4 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS4181127/counts_unfiltered/adata.loom") 
sample4 = sample4[np.isin(sample4.obs.barcode,sample4_barcodes["x"])]
sample4.obs['barcodeID'] = sample4.obs.barcode
for i in range(len(sample4.obs.barcode)):
   sample4.obs['barcodeID'][i] = sample4.obs.barcode[i] + '_7'

sample5_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS4181130.csv")
sample5 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS4181130/counts_unfiltered/adata.loom")
sample5 = sample5[np.isin(sample5.obs.barcode,sample5_barcodes["x"])]
sample5.obs['barcodeID'] = sample5.obs.barcode
for i in range(len(sample5.obs.barcode)):
   sample5.obs['barcodeID'][i] = sample5.obs.barcode[i] + '_8'

sample6_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/LZ007.csv")
sample6 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/LZ007/counts_unfiltered/adata.loom")
sample6 = sample6[np.isin(sample6.obs.barcode,sample6_barcodes["x"])]
sample6.obs['barcodeID'] = sample6.obs.barcode
for i in range(len(sample6.obs.barcode)):
   sample6.obs['barcodeID'][i] = sample6.obs.barcode[i] + '_11'

sample9_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS5883825.csv")
sample9 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS5883825/counts_unfiltered/adata.loom")
sample9 = sample9[np.isin(sample9.obs.barcode,sample6_barcodes["x"])]
sample9.obs['barcodeID'] = sample9.obs.barcode
for i in range(len(sample9.obs.barcode)):
   sample9.obs['barcodeID'][i] = sample9.obs.barcode[i] + '_19'
   
sample10_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS2823408.csv")
sample10 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS2823408/counts_unfiltered/adata.loom")
sample10 = sample10[np.isin(sample10.obs.barcode,sample6_barcodes["x"])]
sample10.obs['barcodeID'] = sample10.obs.barcode
for i in range(len(sample10.obs.barcode)):
   sample10.obs['barcodeID'][i] = sample10.obs.barcode[i] + '_33'

sample11_barcodes = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.barcodes_per_sample/SRS2823409.csv")
sample11 = anndata.read_loom("/project/GorielyLab2021/sbush/ssc_atlas/kb_velo/Homo_sapiens/SRS2823409/counts_unfiltered/adata.loom")
sample11 = sample11[np.isin(sample11.obs.barcode,sample6_barcodes["x"])]
sample11.obs['barcodeID'] = sample11.obs.barcode
for i in range(len(sample11.obs.barcode)):
   sample11.obs['barcodeID'][i] = sample11.obs.barcode[i] + '_34'

sample = sample1.concatenate(sample2, sample3, sample4, sample5, sample6, sample9, sample10, sample11)
sample_index = pd.DataFrame(sample.obs.barcodeID)

clusters = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.clusters.res1.2.csv")
clusters = clusters.rename(columns = {'BC':'barcodeID'})
clusters = clusters.rename(columns = {'CLUSTER':'cluster'})
clusters_ordered = sample_index.merge(clusters, on = "barcodeID")
clusters_ordered = clusters_ordered.iloc[:,1:]
sample.obs['clusters'] = clusters_ordered.values

umap = pd.read_csv("/project/GorielyLab2021/sbush/human_adult_SSCs/results/states0to4.cell_embeddings.csv")
umap = umap.rename(columns = {'Unnamed: 0':'barcodeID'})
umap_ordered = sample_index.merge(umap, on = "barcodeID")
umap_ordered = umap_ordered.iloc[:,1:]
sample.obsm['X_umap'] = umap_ordered.values

scv.pp.filter_and_normalize(sample)
scv.pp.moments(sample)
scv.tl.velocity(sample, mode = "stochastic")
scv.tl.velocity_graph(sample)

scv.pl.velocity_embedding_stream(sample, basis = 'umap', dpi = 300, color = 'clusters', palette={"state 0": "orange", "state 0A/1": "green", "state 0/0B": "red", "state 0B/1": "darkred", "SPG": "yellow", "diff SPG": "pink", "leptotene": "cornflowerblue", "zygotene": "blueviolet"}, save = 'states0to4.velocity_embedding_stream.res1.2.png')
scv.pl.velocity_embedding(sample, basis = 'umap', dpi = 300, color = 'clusters', palette={"state 0": "orange", "state 0A/1": "green", "state 0/0B": "red", "state 0B/1": "darkred", "SPG": "yellow", "diff SPG": "pink", "leptotene": "cornflowerblue", "zygotene": "blueviolet"}, save = 'states0to4.velocity_embedding.res1.2.png')