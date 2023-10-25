# ssc_atlas

This repository contains the scripts used for the analyses described in “Cross-species transcriptomic profiling of spermatogonial stem cells”.

Collectively, they integrate 34 publicly-available post-pubescent human testes scRNA-seq datasets (from [Guo 2018](https://pubmed.ncbi.nlm.nih.gov/30315278/), [Hermann 2018](https://pubmed.ncbi.nlm.nih.gov/30404016/), [Sohni 2019](https://pubmed.ncbi.nlm.nih.gov/30726734/), [Chen 2020](https://doi.org/10.1016/j.gendis.2020.09.004), [Guo 2020](http://pubmed.ncbi.nlm.nih.gov/31928944/), [Shami 2020](https://pubmed.ncbi.nlm.nih.gov/32504559/), [Zhao 2020](https://pubmed.ncbi.nlm.nih.gov/33173058/), [Di Persio 2021](https://pubmed.ncbi.nlm.nih.gov/34622232/), and [Nie 2022](https://pubmed.ncbi.nlm.nih.gov/35504286)) for the purpose of generating a comprehensive transcriptional atlas of human spermatogonial stem cells (hSSCs).

Subsequent scripts parse this atlas to obtain insight into the developmental trajectories of hSSCs by performing comparative analyses of the SSC transcriptome both across human development (using 28 foetal and pre-pubertal samples from ) and across species (using data from [sheep](https://pubmed.ncbi.nlm.nih.gov/33197070/), [pig](https://pubmed.ncbi.nlm.nih.gov/34872612/), [buffalo](https://pubmed.ncbi.nlm.nih.gov/36582818/), [cynomologus macaque](https://pubmed.ncbi.nlm.nih.gov/32795394/), [rat](https://pubmed.ncbi.nlm.nih.gov/35536782/) and mouse (from [Green 2018](https://pubmed.ncbi.nlm.nih.gov/30146481/), [Ernst 2019](https://pubmed.ncbi.nlm.nih.gov/30890697/), [Grive 2019](https://pubmed.ncbi.nlm.nih.gov/30893341), and [Jung 2019](https://pubmed.ncbi.nlm.nih.gov/31237565/))).

The scripts, which should be run in numbered order, perform the following steps and have the following software prerequisites, which should be accessible on the command line. More specific detail on prerequisites is available as comments within the scripts themselves.
