#!/bin/bash

## activate the Kallisto/Bustools environment

conda activate kbpython

## create KB indexes suitable for expression quantification across a range of species

### HUMAN ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Homo_sapiens && cd $_

wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.104.gtf

rm Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.104.gtf

## create a KB index suitable for RNA velocity analysis. See https://www.kallistobus.tools/tutorials/kb_velocity_index/python/kb_velocity_index/

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Homo_sapiens/velocity && cd $_

wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 cdna.fa -f2 intron.fa -c1 cdna_t2c.txt -c2 intron_t2c.txt --workflow lamanno Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.104.gtf

rm Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.104.gtf

### MACAQUE (CYNOMOLGUS) ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Macaca_fascicularis && cd $_

wget http://ftp.ensembl.org/pub/release-107/gtf/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/macaca_fascicularis/cdna/Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf

rm Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa Macaca_fascicularis.Macaca_fascicularis_6.0.dna.toplevel.fa Macaca_fascicularis.Macaca_fascicularis_6.0.107.gtf

### MACAQUE (RHESUS) ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Macaca_mulatta && cd $_

wget http://ftp.ensembl.org/pub/release-107/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.107.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Macaca_mulatta.Mmul_10.cdna.all.fa Macaca_mulatta.Mmul_10.dna.toplevel.fa Macaca_mulatta.Mmul_10.107.gtf

rm Macaca_mulatta.Mmul_10.cdna.all.fa Macaca_mulatta.Mmul_10.dna.toplevel.fa Macaca_mulatta.Mmul_10.107.gtf

### MOUSE ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Mus_musculus && cd $_

wget http://ftp.ensembl.org/pub/release-107/gtf/mus_musculus/Mus_musculus.GRCm39.107.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Mus_musculus.GRCm39.cdna.all.fa Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus.GRCm39.107.gtf

rm Mus_musculus.GRCm39.cdna.all.fa Mus_musculus.GRCm39.dna.primary_assembly.fa Mus_musculus.GRCm39.107.gtf

### PIG ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Sus_scrofa && cd $_

wget http://ftp.ensembl.org/pub/release-107/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.107.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz
gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Sus_scrofa.Sscrofa11.1.cdna.all.fa Sus_scrofa.Sscrofa11.1.dna.toplevel.fa Sus_scrofa.Sscrofa11.1.107.gtf

rm Sus_scrofa.Sscrofa11.1.cdna.all.fa Sus_scrofa.Sscrofa11.1.dna.toplevel.fa Sus_scrofa.Sscrofa11.1.107.gtf

### RAT ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Rattus_norvegicus && cd $_

wget https://ftp.ensembl.org/pub/release-108/gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.108.gtf.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-108/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.mRatBN7.2.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Rattus_norvegicus.mRatBN7.2.cdna.all.fa Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa Rattus_norvegicus.mRatBN7.2.108.gtf

rm Rattus_norvegicus.mRatBN7.2.108.gtf Rattus_norvegicus.mRatBN7.2.cdna.all.fa Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa

### SHEEP ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Ovis_aries && cd $_

wget http://ftp.ensembl.org/pub/release-107/gtf/ovis_aries_rambouillet/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.107.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/ovis_aries_rambouillet/dna/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/ovis_aries_rambouillet/cdna/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa Ovis_aries_rambouillet.Oar_rambouillet_v1.0.107.gtf

rm Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa Ovis_aries_rambouillet.Oar_rambouillet_v1.0.107.gtf

### BUFFALO ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Bubalus_bubalis && cd $_

wget https://ftp.ensembl.org/pub/rapid-release/species/Bubalus_bubalis/GCA_003121395.1/ensembl/geneset/2020_06/Bubalus_bubalis-GCA_003121395.1-2020_06-genes.gtf.gz
wget https://ftp.ensembl.org/pub/rapid-release/species/Bubalus_bubalis/GCA_003121395.1/ensembl/genome/Bubalus_bubalis-GCA_003121395.1-unmasked.fa.gz
wget https://ftp.ensembl.org/pub/rapid-release/species/Bubalus_bubalis/GCA_003121395.1/ensembl/geneset/2020_06/Bubalus_bubalis-GCA_003121395.1-2020_06-cdna.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Bubalus_bubalis-GCA_003121395.1-2020_06-cdna.fa Bubalus_bubalis-GCA_003121395.1-unmasked.fa Bubalus_bubalis-GCA_003121395.1-2020_06-genes.gtf

rm Bubalus_bubalis-GCA_003121395.1-2020_06-cdna.fa Bubalus_bubalis-GCA_003121395.1-unmasked.fa Bubalus_bubalis-GCA_003121395.1-2020_06-genes.gtf

### PANDA ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Ailuropoda_melanoleuca && cd $_

wget https://ftp.ensembl.org/pub/release-109/gtf/ailuropoda_melanoleuca/Ailuropoda_melanoleuca.ASM200744v2.109.gtf.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/ailuropoda_melanoleuca/dna/Ailuropoda_melanoleuca.ASM200744v2.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/ailuropoda_melanoleuca/cdna/Ailuropoda_melanoleuca.ASM200744v2.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Ailuropoda_melanoleuca.ASM200744v2.cdna.all.fa Ailuropoda_melanoleuca.ASM200744v2.dna.toplevel.fa Ailuropoda_melanoleuca.ASM200744v2.109.gtf

rm Ailuropoda_melanoleuca.ASM200744v2.cdna.all.fa Ailuropoda_melanoleuca.ASM200744v2.dna.toplevel.fa Ailuropoda_melanoleuca.ASM200744v2.109.gtf

### CHICKEN ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Gallus_gallus && cd $_

wget https://ftp.ensembl.org/pub/release-109/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/cdna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf

rm Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf

### BABOON ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Papio_anubis && cd $_

wget https://ftp.ensembl.org/pub/release-109/gtf/papio_anubis/Papio_anubis.Panubis1.0.109.gtf.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/papio_anubis/dna/Papio_anubis.Panubis1.0.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/papio_anubis/cds/Papio_anubis.Panubis1.0.cds.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Papio_anubis.Panubis1.0.cds.all.fa Papio_anubis.Panubis1.0.dna.toplevel.fa Papio_anubis.Panubis1.0.109.gtf

rm Papio_anubis.Panubis1.0.cds.all.fa Papio_anubis.Panubis1.0.dna.toplevel.fa Papio_anubis.Panubis1.0.109.gtf

### YAK ###

mkdir /project/GorielyLab2021/sbush/ssc_atlas/indexes/Bos_grunniens && cd $_

wget https://ftp.ensembl.org/pub/release-109/gtf/bos_grunniens/Bos_grunniens.LU_Bosgru_v3.0.109.gtf.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/bos_grunniens/dna/Bos_grunniens.LU_Bosgru_v3.0.dna.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/release-109/fasta/bos_grunniens/cds/Bos_grunniens.LU_Bosgru_v3.0.cds.all.fa.gz

gunzip *.gz

kb ref -i index.idx -g t2g.txt -f1 Bos_grunniens.LU_Bosgru_v3.0.cds.all.fa Bos_grunniens.LU_Bosgru_v3.0.dna.toplevel.fa Bos_grunniens.LU_Bosgru_v3.0.109.gtf

rm Bos_grunniens.LU_Bosgru_v3.0.cds.all.fa Bos_grunniens.LU_Bosgru_v3.0.dna.toplevel.fa Bos_grunniens.LU_Bosgru_v3.0.109.gtf