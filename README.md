# spg_atlas

This repository contains the scripts used for the analyses described in Bush _et al._ (2024) [Adult Human, but Not Rodent, Spermatogonial Stem Cells Retain States with a Foetal-like Signature](https://pubmed.ncbi.nlm.nih.gov/38727278/).

Collectively, they integrate 34 publicly-available post-pubescent human testes scRNA-seq datasets (from [Guo 2018](https://pubmed.ncbi.nlm.nih.gov/30315278/), [Hermann 2018](https://pubmed.ncbi.nlm.nih.gov/30404016/), [Sohni 2019](https://pubmed.ncbi.nlm.nih.gov/30726734/), [Chen 2020](https://doi.org/10.1016/j.gendis.2020.09.004), [Guo 2020](http://pubmed.ncbi.nlm.nih.gov/31928944/), [Shami 2020](https://pubmed.ncbi.nlm.nih.gov/32504559/), [Zhao 2020](https://pubmed.ncbi.nlm.nih.gov/33173058/), [Di Persio 2021](https://pubmed.ncbi.nlm.nih.gov/34622232/), and [Nie 2022](https://pubmed.ncbi.nlm.nih.gov/35504286)) for the purpose of generating a comprehensive transcriptional atlas of human spermatogonia, accessible via ShinyApp [here](https://trainingidn.shinyapps.io/hu_spermatogonial_atlas_shinyApp/). Subsequent scripts parse this atlas to obtain insight into the developmental trajectories of spermatogonia (SPG) by performing comparative analyses of the SPG transcriptome both across human development (using 28 foetal and pre-pubertal samples from [Guo 2018](https://pubmed.ncbi.nlm.nih.gov/30315278/), [Sohni 2019](https://pubmed.ncbi.nlm.nih.gov/30726734/), [Guo 2020](http://pubmed.ncbi.nlm.nih.gov/31928944/), [Zhao 2020](https://pubmed.ncbi.nlm.nih.gov/33173058/), [Guo 2021](https://pubmed.ncbi.nlm.nih.gov/33453151/), [Voigt 2022](https://pubmed.ncbi.nlm.nih.gov/35856882/), and [Wang 2022](https://pubmed.ncbi.nlm.nih.gov/35513251/)) and across species (using data from [sheep](https://pubmed.ncbi.nlm.nih.gov/33197070/), [pig](https://pubmed.ncbi.nlm.nih.gov/34872612/), [buffalo](https://pubmed.ncbi.nlm.nih.gov/36582818/), [cynomologus macaque](https://pubmed.ncbi.nlm.nih.gov/32795394/), [rat](https://pubmed.ncbi.nlm.nih.gov/35536782/) and mouse (from [Green 2018](https://pubmed.ncbi.nlm.nih.gov/30146481/), [Ernst 2019](https://pubmed.ncbi.nlm.nih.gov/30890697/), [Grive 2019](https://pubmed.ncbi.nlm.nih.gov/30893341), and [Jung 2019](https://pubmed.ncbi.nlm.nih.gov/31237565/))).

The scripts, which should be run in numbered order, perform the following steps and have the following software prerequisites, which should be accessible on the command line. More specific detail on prerequisites is available as comments within the scripts themselves.

## 1.download_fqs.pl

This script contains commands for downloading raw sequencing data from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home), either by obtaining fqs directly with [Aspera](https://www.ibm.com/products/aspera) or, if raw 10X data was instead uploaded in BAM format, converting it using [bamtofastq](https://support.10xgenomics.com/docs/bamtofastq).

## 2.make_kb_index.sh

This script contains commands for creating two [Kallisto/Bustools](https://www.kallistobus.tools/) (KB) indexes, one for expression quantification and one for RNA velocity analysis. These commands are as detailed in the [KB tutorial](https://www.kallistobus.tools/tutorials/kb_velocity_index/python/kb_velocity_index/). The output is an index file, index.idx, as well as a transcript-to-gene lookup table, t2g.txt.

## 3.update_t2g_to_account_for_duplicate_gene_ids.pl

This script revises the t2g.txt file output by the KB index-building command. We will ultimately be using KB to quantify expression on a per-transcript basis, and then aggregating results to gene level. For this purpose, we wish to ensure each Ensembl transcript ID maps to a unique gene name. To do this, we determine whether a given gene name is assigned to multiple Ensembl gene IDs, revising those that do to the form “gene name/Ensembl gene ID”, e.g. “Y_RNA/ENSG00000206677”, “Y_RNA/ENSG00000206646” and “Y_RNA/ENSG00000202272.”

## 4.run_kb.pl

This script contains commands for running the KB ‘count’ workflow to quantify expression per sample, producing (among other output) a cell x gene count matrix (as an RDS object) and a [UMAP](https://ui.adsabs.harvard.edu/abs/2018arXiv180203426M) visualisation.

To do this, the script calls [R](https://www.r-project.org/) to perform a number of QC steps as well as to process each count matrix using [Seurat](https://satijalab.org/seurat/), implementing commands adapted from those given in the [KB atlas-building user guide](https://www.kallistobus.tools/tutorials/kb_building_atlas/r/kb_analysis_0_r/).

Note that as an initial QC step for each sample, KB automatically filters cells on the basis of a [knee plot](https://www.cell.com/fulltext/S0092-8674(15)00549-8), retaining only those above a sample-specific, empirically-derived, inflection point. Parameters hard-coded into this script also retain only those genes detected in at least 3 cells, and only those cells where the total number of genes was > 1000 and < 10,000, the total number of UMIs was > 2000 and < 50,000 and the proportion of mitochondrial genes detected was < 5% of the total.

The Seurat commands follow developer recommendations and implement the normalisation procedure [SCTransform](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1), setting method to ‘glmGamPoi’ and selecting mitochondrial gene content, total gene count, total RNA count, and both ‘S score’ and ‘G2M score’ as regression variables. Seurat’s FindClusters function is also configured to run across an intentionally large range of resolutions, from 0 to 2.0 at intervals of 0.2. The resolution parameter sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. [Seurat guidelines](https://satijalab.org/seurat/archive/v3.1/pbmc3k_tutorial.html) on this issue are that a broad range of values are generally acceptable (producing ‘typically good’ results) but that optimisation is certainly possible:

>we find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.

To aid in optimisation, the script produces a [Clustree](https://academic.oup.com/gigascience/article/7/7/giy083/5052205) dendrogram, showing the impact of different resolution values on the number and composition of the clusters.

At increasingly high resolutions, we expect to see this dendrogram becoming messier, comprising nodes with multiple incoming edges. This indicates over-clustering: smaller-sized clusters becoming increasingly fragmented as a proportion of their cells are partitioned into ever-smaller, separate, clusters. Conversely, if the resolution parameter is too low, the data will lack specificity. A more detailed explanation of the problem is available [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-03957-4):

>without foreknowledge of cell types, it is hard to address the quality of the chosen clusters, and whether the cells have been under- or over-clustered. In general, under-clustering occurs when clusters are too broad and mask underlying biological structure. Near-optimal clustering is when most clusters relate to known or presumed cell types, with relevant biological distinctions revealed and without noisy, unreliable, or artifactual sub-populations. When cells are slightly over-clustered, non-relevant subdivisions have been introduced; however, these subclusters can still be merged to recover appropriate cell types. Once severe over-clustering occurs, however, some clusters may be shattered, meaning they are segregated based on non-biological variation to the point where iterative re-merging cannot recover the appropriate cell types.

## 5.integrate_all_whole_testes_datasets.R

This script ingests the set of RDS objects created by script 4, one per sample, and applies Seurat’s [rPCA](https://satijalab.org/seurat/articles/integration_rpca.html) integration strategy. This is a memory-saving technique [recommended for large datasets](https://satijalab.org/seurat/articles/integration_large_datasets.html), requiring the user define a ‘mutual neighbourhood’ in which to search for integration anchors. For this purpose, we use as our ‘neighbourhood’ the combined set of [Di Persio](https://pubmed.ncbi.nlm.nih.gov/34622232/), [Sohni](https://pubmed.ncbi.nlm.nih.gov/30726734/) and [Zhao](https://pubmed.ncbi.nlm.nih.gov/33173058/) samples, as together these contain the greatest number of cells (for number of cells per sample, see **Supplementary Table 3**; this is produced by script 7 in this repo). This script also produces UMAP plots of the integrated atlas, coloured according to sample ID, source and donor age, and demonstrating that there are no discernible batch effects introduced by either variable (**Supplementary Figure 2**), as well as a Clustree dendrogram (see above).
In the next script we will re-run the clustering having decided on an optimal resolution value (by inspection of both the Clustree plot and the cluster boundaries using cellxgene), alongside cleaning the data object and exporting finalised metadata. Cellxgene takes as input biomarkers.csv, a manually curated list of testis cell type and spermatogonial stem cell state-related genes, and used to inform annotation (for sources, see **Supplementary Table 4**).

## 6.refine_clustering_of_whole_testes_atlas.R

This script ingests the ‘first pass’ integrated dataset created by script 6, to which a range of cluster resolutions had been applied. It now re-clusters the atlas using a single, empirically derived, resolution value (in this case, 0.15), producing a refined UMAP. This script also runs a global differential gene expression analysis (i.e. Seurat’s [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers) function) to facilitate an unbiased assignment of cell type identities. The results of the differential expression analysis are found in **Supplementary Table 5**.

## 7.count_num_cells_and_genes_passing_qc.pl

This script parses the output of scripts 4 and 6 to produce a summary table: the number, and percentage, of cells and genes per sample that pass QC thresholds (e.g. number of genes, number of UMIs, and % mitochondria content per cell) as well as the average number of genes/cell that each sample contributes to the overall atlas. This table is **Supplementary Table 3**.

## 8.subset_whole_testes_atlas_to_obtain_States0to4.R

This script ingests the whole testes atlas (produced by script 6), subsets the undifferentiated and differentiating SPG clusters (i.e. cells near or at the onset of the meiotic program), and then re-integrates all cells to produce a ‘SPG atlas’ (as an RDS object). Note that unlike script 5, the integration strategy applied in this case is [CCA](https://satijalab.org/seurat/articles/integration_rpca.html) (canonical correlation analysis), not rPCA, as this is instead
>well-suited for identifying anchors when cell types are conserved, but there are very substantial differences in gene expression across experiments

and so is considered appropriate for this dataset: the majority of cells are of the same type (SPG), differing only in transcriptional state, and show substantial differences in gene expression across the dataset (because the dataset represents differentiation and commitment to meiosis).
While the whole-testes atlas was created using two scripts (numbers 5 and 6), this script implements the same ‘two-pass’ approach, but in one place.

Note that the SPG data were initially clustered at very high resolution (5.0) such that cluster IDs could be assigned to three small contaminating cell populations; these are removed and the dataset iteratively re-integrated to converge on the final atlas. The final SPG atlas is illustrated in **Figure 1**, with the undifferentiated SPGs forming a ‘ring’ structure (bottom right of figure) from which the differentiation trajectory progresses upwards to the onset of the meiotic program (top left of figure). This ‘ring’ can also be seen in [Figure 2C of Sohni 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6402825/figure/F2/), [Figure 5C of Di Persio 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8484693/figure/fig5/), and [Figure S1 of Huang 2023](https://static-content.springer.com/esm/art%3A10.1038%2Fs41556-023-01232-7/MediaObjects/41556_2023_1232_MOESM1_ESM.pdf), although in each case its biological significance is unclear. To gain insight into the transcriptional states comprising this ring, this script clusters the SPG atlas at three different resolutions, empirically selected as the minimum values necessary to partition the ring into two, three or four compartments (resolutions 0.3, 1.1 and 1.2, respectively, illustrated in **Supplementary Figure 9**).

## 9.run_diff_exp_analysis_on_States0to4.R

This script contains commands for running Seurat’s [FindAllMarkers](https://satijalab.org/seurat/reference/findallmarkers) and [FindMarkers](https://satijalab.org/seurat/reference/findmarkers) functions, the purpose of which are to identify putative discriminative biomarkers of each cluster, either relative to all other clusters or in a pairwise comparison to one specified other cluster, respectively.

## 10.run_SPIA_on_States0to4.pl

This script parses the differential expression analysis results (from script 9) to create the input files necessary to run the KEGG enrichment analysis package [SPIA](https://bioconductor.org/packages/SPIA/). R code for running this package is given as comments within the script. The results of the analysis are given as **Supplementary Table 12**.

## 11.count_num_of_cells_per_cluster_in_which_each_gene_is_expressed.R

This script parses the RDS object created in script 8 and outputs, per cluster and per Seurat resolution, tables showing the number of cells in that cluster expressing each gene, the proportion of cells in that cluster expressing each gene, and the average expression per gene of all cells in which that gene is detectably expressed. These tables are used as input for the next script, which aggregates this data into one large ‘atlas’ table.

## 12.create_summary_table_of_atlas

_12a.create_summary_table_of_whole_testes_atlas.pl_

_12b.create_summary_table_of_States0to4_resolution_0.3.pl_

_12c.create_summary_table_of_States0to4_resolution_1.1.pl_

_12d.create_summary_table_of_States0to4_resolution_1.2.pl_

These scripts aggregate the various atlas metadata files to create summary tables for both the whole testis and SPG atlases (the latter at three different resolutions) which detail, per gene, its average expression across all cells per cell type/SPG state, the proportion of cells in which it is detectably expressed, and Seurat’s all-against-all differential expression analysis results (from script 9). The intention with these tables is to facilitate rapid identification of prospective biomarkers for a given cluster (i.e. cell type/state). The tables are provided as **Supplementary Tables 5, 6, 7 and 8**.

## 13.turn_go_terms_into_map_file.pl

This script parses a per-gene list of GO term accessions downloaded from [Ensembl BioMart](https://www.ensembl.org/info/data/biomart/index.html), to produce a ‘.map’ file suitable for use with the R package [topGO](https://bioconductor.org/packages/topGO/).

## 14.run_topgo_on_cluster_markers_for_atlas

_14a.run_topgo_on_cluster_markers_for_whole_testes_atlas.pl_

_14b.run_topgo_on_cluster_markers_for_SPG_atlas.pl_

_14c.run_topgo_on_cluster_markers_for_SPG_atlas_in_X_vs_Y_comparison.pl_

These scripts parse both the differential expression analysis results (from script 9) and the per-gene GO term list (from script 13) to create the input files and R code necessary to run the GO term enrichment analysis package topGO.

## 15.parse_topgo_results_for_atlas

_15a.parse_topgo_results_for_whole_testes_atlas.pl_

_15b.parse_topgo_results_for_SPG_atlas.pl_

_15c.parse_topgo_results_for_SPG_atlas_in_X_vs_Y_comparison.pl_

These scripts parse the topGO output into user-friendly results tables, alongside providing additional summary statistics (e.g. the percentage of the total genes with each GO term present in each cluster). The GO terms enriched among the set of genes differentially expressed in each cluster in each atlas, for both all-against-all and pairwise comparisons, are given as **Supplementary Tables 9 to 11**.

## 16.run_scVelo

_16a.parse_States0to4_to_obtain_only_reads_from_particular_cells.pl_

_16b.run_scvelo_on_States0to4_atlas_at_resolution1_1.py_

_16c.run_scvelo_on_States0to4_atlas_at_resolution1_2.py_

To calculate RNA velocity, we first need to obtain all reads associated with the cell barcodes used in the SPG atlas. We’ll only need to calculate velocity for these cells, because it’s only those that will later be projected onto the atlas UMAP. To do this, script 16a exports the list of cell barcodes and their associated sample ID from the atlas RDS, created by script 8. This is later used as input for scripts 16b and 16c, which contain the [scVelo](https://scvelo.readthedocs.io/) commands for creating the velocity plots of **Figure 2** and **Supplementary Figure 11**.

## 17.create_SpeciesX_to_human_symbol_lookup_table.pl

This script parses one-to-one gene orthology data downloaded from [Ensembl BioMart](https://www.ensembl.org/info/data/biomart/index.html) to create a series of lookup tables, one per species, showing the equivalent human name for each gene. These tables are later used when projecting data from each species onto the human SPG UMAP.

## 18.project_data_onto_atlas

_18a.project_mouse_time_course_onto_adult_whole_testes.R_

_18b.project_mouse_germline_time_course_onto_States0to4_resolution1_1.R_

_18c.project_human_time_course_onto_adult_whole_testes.R_

_18d.project_human_germline_time_course_onto_States0to4_resolution1_1.R_

_18e.project_sheep_pig_buffalo_and_cynomolgus_macaque_onto_adult_whole_testes.R_

_18f.project_sheep_pig_buffalo_and_cynomolgus_macaque_onto_States0to4_resolution1_1.R_

_18g.project_rat_onto_adult_whole_testes.R_

_18h.project_rat_onto_States0to4_resolution1_1.R_

_18i.project_Hermann_samples_onto_States0to4_resolution1_1_to_validate_projection_strategy.R_

_18j.project_undiff_mouse_SSCs_onto_States0to4_resolution1_1.R_

These scripts contain the Seurat commands necessary to project cells from a given sample onto the SPG atlas. They employ a gating strategy to first distinguish somatic from germline cells, and then re-project only the latter onto the SPG atlas UMAP. These scripts create **Figure 3** (sheep, pig, buffalo, macaque), **Figure 4** (mouse), and **Figure 5** (rat).

## 19.count_num_of_cells_projected_to_each_State0to4_cluster_for_human/mouse/other

_19a.count_num_of_cells_projected_to_each_State0to4_cluster_for_human_germline_time_course.pl_

_19b.count_num_of_cells_projected_to_each_State0to4_cluster_for_mouse_time_course.pl_

_19c.count_num_of_cells_projected_to_each_State0to4_cluster_for_all_other_species.pl_

These scripts parse the output of the ‘projection’ scripts to count the proportion of germ cells at each time point classified as state 0, 0A/1, 0B, and so on, as shown in **Supplementary Tables 15, 16, and 17**. These tables constitute the raw data from which **Figures 2A** (human) and **2B** (mouse) are created.

## 20.parse_summary_table_of_States0to4_to_create_shortlist_of_quiescence_associated_genes.pl

This script parses the atlas summary table (from script 12) to create a shortlist of quiescence-associated genes and their relative expression in each of the undifferentiated SPG clusters; this is **Supplementary Table 14**.
