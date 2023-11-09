=head
BEFORE USAGE:

# the commands in this script automate the download of fqs/bams from the ENA, where possible
# note that for some species and samples, we must download data manually

## for the HUMAN data, we must download one SSC sample manually

# data from Chen 2020 is not mirrored on the ENA and should be manually downloaded instead from the NODE database (National Omics Data Encyclopaedia), as follows
# this study comprises 2 samples, one OA (https://www.biosino.org/node/sample/detail/OES026004 | https://www.biosino.org/node/run/detail/OER039859) and one NOA (https://www.biosino.org/node/sample/detail/OES026003 | https://www.biosino.org/node/run/detail/OER039858)

cd /project/GorielyLab2021/sbush/ssc_atlas/fq/Homo_sapiens/
mkdir OES026004 && cd $_
wget https://www.biosino.org/download/node/data/public/OED111526 --tries 0 -O OES026004.1.fq.gz
wget https://www.biosino.org/download/node/data/public/OED113871 --tries 0 -O OES026004.2.fq.gz
mkdir OES026003 && cd $_
wget https://www.biosino.org/download/node/data/public/OED111786 --tries 0 -O OES026003.1.fq.gz
wget https://www.biosino.org/download/node/data/public/OED111785 --tries 0 -O OES026003.2.fq.gz

## for the CYNOMOLGUS data, we must download all samples manually

# infant
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/infant
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/infant
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S1_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S1_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S33_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S33_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S34_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S34_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S35_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S35_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S36_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S36_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S41_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S41_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S42_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S42_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S43_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S43_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S44_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Infant_S44_L004_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > infant.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > infant.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# juvenile
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/juvenile
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/juvenile
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S1_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S1_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S3_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S3_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S45_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S45_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S46_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S46_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S47_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S47_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S48_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S48_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S4_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S4_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S5_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S5_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S6_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Juvenile_S6_L004_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > juvenile.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > juvenile.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# adult 1
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/adult1
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/adult1
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S1_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S1_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S29_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S29_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S30_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S30_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S31_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S31_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S32_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S32_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S37_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S37_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S38_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S38_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S39_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S39_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S40_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult1_S40_L004_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > adult1.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > adult1.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# adult 2
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/adult2
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/Macaca_fascicularis/fq/adult2
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S12_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S12_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S13_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S13_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S14_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S14_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S15_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S15_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S1_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S1_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S28_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S28_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S7_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S7_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S8_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S8_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S9_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8979/Adult2_S9_L004_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > adult2.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > adult2.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

## for the MOUSE data, we must download some samples manually

# ERS2575676
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575676
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575676
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S5_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S5_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S5_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S5_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S6_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S6_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S6_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S6_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S7_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S7_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S7_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S7_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S8_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15983_S8_L006_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575676.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575676.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575677
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575677
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575677
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S1_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S1_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S1_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S1_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S2_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S2_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S2_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S2_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S3_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S3_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S3_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S3_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S4_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S4_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S4_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do15984_S4_L007_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575677.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575677.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575678
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575678
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575678
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S1_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S2_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S3_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17622_S4_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575678.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575678.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575679
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575679
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575679
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S5_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S6_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S7_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17623_S8_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575679.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575679.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575680
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575680
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575680
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S17_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S18_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S19_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17815_S20_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575680.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575680.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575681
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575681
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575681
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S13_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S14_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S15_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17816_S16_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575681.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575681.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575682
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575682
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575682
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S10_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S11_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S12_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17821_S9_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575682.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575682.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575683
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575683
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575683
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S29_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S30_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S31_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17824_S32_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575683.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575683.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575684
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575684
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575684
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S21_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S22_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S23_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17825_S24_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575684.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575684.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575685
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575685
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575685
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S29_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S30_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S31_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do17827_S32_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575685.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575685.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575686
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575686
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575686
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S17_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S18_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S19_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18195_S20_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575686.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575686.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575687
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575687
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575687
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S13_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S14_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S15_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18196_S16_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575687.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575687.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575688
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575688
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575688
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S25_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S26_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S27_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18197_S28_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575688.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575688.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575689
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575689
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575689
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S21_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S22_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S23_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18198_S24_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575689.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575689.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS2575690
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575690
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS2575690
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S29_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S30_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S31_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L001_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L001_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L002_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L002_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L003_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L003_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L004_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L004_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L005_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L005_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L006_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L006_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L007_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L007_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do18199_S32_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS2575690.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS2575690.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS3000379
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS3000379
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS3000379
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S1_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S1_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S2_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S2_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S3_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S3_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S4_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26386_S4_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS3000379.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS3000379.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

# ERS3000380
mkdir /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS3000380
cd /t1-data/project/GorielyLab2021/sbush/ssc_atlas/fq/Mus_musculus/ERS3000380
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S5_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S5_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S6_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S6_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S7_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S7_L008_R2_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S8_L008_R1_001.fastq.gz
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6946/do26387_S8_L008_R2_001.fastq.gz
zcat *_R1_001.fastq.gz | gzip > ERS3000380.1.fq.gz
zcat *_R2_001.fastq.gz | gzip > ERS3000380.2.fq.gz
rm *_R1_001.fastq.gz
rm *_R2_001.fastq.gz

## for the PANDA data, we must download all samples manually (because the authors have split up their paired-end files and uploaded them to the ENA separately)

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Ailuropoda_melanoleuca/SRS9879449 && cd $_
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/078/SRR15569478/SRR15569478.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/080/SRR15569480/SRR15569480.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/082/SRR15569482/SRR15569482.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/084/SRR15569484/SRR15569484.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/088/SRR15569488/SRR15569488.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/094/SRR15569494/SRR15569494.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/068/SRR15569468/SRR15569468.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/073/SRR15569473/SRR15569473.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/079/SRR15569479/SRR15569479.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/087/SRR15569487/SRR15569487.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/090/SRR15569490/SRR15569490.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/092/SRR15569492/SRR15569492.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/093/SRR15569493/SRR15569493.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/095/SRR15569495/SRR15569495.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/096/SRR15569496/SRR15569496.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/098/SRR15569498/SRR15569498.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/067/SRR15569467/SRR15569467.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/069/SRR15569469/SRR15569469.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/070/SRR15569470/SRR15569470.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/071/SRR15569471/SRR15569471.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/072/SRR15569472/SRR15569472.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/074/SRR15569474/SRR15569474.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/075/SRR15569475/SRR15569475.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/076/SRR15569476/SRR15569476.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/077/SRR15569477/SRR15569477.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/081/SRR15569481/SRR15569481.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/083/SRR15569483/SRR15569483.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/085/SRR15569485/SRR15569485.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/086/SRR15569486/SRR15569486.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/089/SRR15569489/SRR15569489.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/091/SRR15569491/SRR15569491.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/097/SRR15569497/SRR15569497.fastq.gz

zcat SRR15569468.fastq.gz SRR15569470.fastq.gz SRR15569472.fastq.gz SRR15569474.fastq.gz SRR15569486.fastq.gz SRR15569477.fastq.gz SRR15569479.fastq.gz SRR15569481.fastq.gz SRR15569483.fastq.gz SRR15569485.fastq.gz SRR15569488.fastq.gz SRR15569490.fastq.gz SRR15569492.fastq.gz SRR15569494.fastq.gz SRR15569496.fastq.gz SRR15569498.fastq.gz | gzip > SRS9879449.1.fq.gz
zcat SRR15569467.fastq.gz SRR15569469.fastq.gz SRR15569471.fastq.gz SRR15569473.fastq.gz SRR15569475.fastq.gz SRR15569476.fastq.gz SRR15569478.fastq.gz SRR15569480.fastq.gz SRR15569482.fastq.gz SRR15569484.fastq.gz SRR15569487.fastq.gz SRR15569489.fastq.gz SRR15569491.fastq.gz SRR15569493.fastq.gz SRR15569495.fastq.gz SRR15569497.fastq.gz | gzip > SRS9879449.2.fq.gz

## for the CHICKEN data, we must download all samples manually

# we can see at https://www.ebi.ac.uk/ena/browser/view/PRJNA761874 that although the data are paired-end, instead of the expected two fastq files, there is only one. From the uploaded file name ('run alias') it appears that this is the index file (I1_001.fastq.gz) rather than R1 or R2.
# so where are the R1 and R2 files? They appear to have been lost in the transition from 'original format' data to 'SRA archive data'
# if we look up the 'original format' files we can see that the R1 and R2 do, in fact, exist: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR15827453&display=data-access
# the way to resolve this problem is to download via the SRA Toolkit instead of directly

module add sratoolkit/3.0.0

# E2.5
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099509 && cd $_
prefetch SRR15827463
fastq-dump --split-files --gzip SRR15827463/SRR15827463.sra
rm -r SRR15827463
mv SRR15827463_2.fastq.gz SRS10099509.1.fq.gz
mv SRR15827463_3.fastq.gz SRS10099509.2.fq.gz
rm SRR15827463_1.fastq.gz

# E6
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099510 && cd $_
prefetch SRR15827462
fastq-dump --split-files --gzip SRR15827462/SRR15827462.sra
rm -r SRR15827462
mv SRR15827462_2.fastq.gz SRS10099510.1.fq.gz
mv SRR15827462_3.fastq.gz SRS10099510.2.fq.gz
rm SRR15827462_1.fastq.gz

# E8
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099514 && cd $_
prefetch SRR15827458
fastq-dump --split-files --gzip SRR15827458/SRR15827458.sra
rm -r SRR15827458
mv SRR15827458_2.fastq.gz SRS10099514.1.fq.gz
mv SRR15827458_3.fastq.gz SRS10099514.2.fq.gz
rm SRR15827458_1.fastq.gz

# E12
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099515 && cd $_
prefetch SRR15827457
fastq-dump --split-files --gzip SRR15827457/SRR15827457.sra
rm -r SRR15827457
mv SRR15827457_2.fastq.gz SRS10099515.1.fq.gz
mv SRR15827457_3.fastq.gz SRS10099515.2.fq.gz
rm SRR15827457_1.fastq.gz

# E16
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099516 && cd $_
prefetch SRR15827456
fastq-dump --split-files --gzip SRR15827456/SRR15827456.sra
rm -r SRR15827456
mv SRR15827456_2.fastq.gz SRS10099516.1.fq.gz
mv SRR15827456_3.fastq.gz SRS10099516.2.fq.gz
rm SRR15827456_1.fastq.gz

# hatch
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099517 && cd $_
prefetch SRR15827455
fastq-dump --split-files --gzip SRR15827455/SRR15827455.sra
rm -r SRR15827455
mv SRR15827455_2.fastq.gz SRS10099517.1.fq.gz
mv SRR15827455_3.fastq.gz SRS10099517.2.fq.gz
rm SRR15827455_1.fastq.gz

# 1 week post-hatching
mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Gallus_gallus/SRS10099518 && cd $_
prefetch SRR15827454
fastq-dump --split-files --gzip SRR15827454/SRR15827454.sra
rm -r SRR15827454
mv SRR15827454_2.fastq.gz SRS10099518.1.fq.gz
mv SRR15827454_3.fastq.gz SRS10099518.2.fq.gz
rm SRR15827454_1.fastq.gz

## for the Singh 2023 MACAQUE data, we must download all samples manually

module add sratoolkit/3.0.0

# Rhesus37136A

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Macaca_mulatta/SRS16347199 && cd $_
prefetch --max-size 100G SRR22963910
fastq-dump --split-files --gzip SRR22963910/SRR22963910.sra
rm -r SRR22963910
mv SRR22963910_2.fastq.gz SRS16347199.1.fq.gz
mv SRR22963910_3.fastq.gz SRS16347199.2.fq.gz
rm SRR22963910_1.fastq.gz

# Rhesus37136B

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Macaca_mulatta/SRS16347200 && cd $_
prefetch --max-size 100G SRR22963909
fastq-dump --split-files --gzip SRR22963909/SRR22963909.sra
rm -r SRR22963909
mv SRR22963909_2.fastq.gz SRS16347200.1.fq.gz
mv SRR22963909_3.fastq.gz SRS16347200.2.fq.gz
rm SRR22963909_1.fastq.gz

# Rhesus34403

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Macaca_mulatta/SRS16347195 && cd $_
prefetch --max-size 100G SRR22963913 SRR22963914 SRR22963916 SRR22963917 SRR22963918 SRR22963919 SRR22963920 SRR22963921 SRR22963922 SRR22963923 SRR22963924 SRR22963925 SRR22963926 SRR22963927 SRR22963928
fastq-dump --split-files --gzip SRR22963913/SRR22963913.sra
fastq-dump --split-files --gzip SRR22963914/SRR22963914.sra
fastq-dump --split-files --gzip SRR22963916/SRR22963916.sra
fastq-dump --split-files --gzip SRR22963917/SRR22963917.sra
fastq-dump --split-files --gzip SRR22963918/SRR22963918.sra
fastq-dump --split-files --gzip SRR22963919/SRR22963919.sra
fastq-dump --split-files --gzip SRR22963920/SRR22963920.sra
fastq-dump --split-files --gzip SRR22963921/SRR22963921.sra
fastq-dump --split-files --gzip SRR22963922/SRR22963922.sra
fastq-dump --split-files --gzip SRR22963923/SRR22963923.sra
fastq-dump --split-files --gzip SRR22963924/SRR22963924.sra
fastq-dump --split-files --gzip SRR22963925/SRR22963925.sra
fastq-dump --split-files --gzip SRR22963926/SRR22963926.sra
fastq-dump --split-files --gzip SRR22963927/SRR22963927.sra
fastq-dump --split-files --gzip SRR22963928/SRR22963928.sra
rm -r SRR22963913
rm -r SRR22963914
rm -r SRR22963916
rm -r SRR22963917
rm -r SRR22963918
rm -r SRR22963919
rm -r SRR22963920
rm -r SRR22963921
rm -r SRR22963922
rm -r SRR22963923
rm -r SRR22963924
rm -r SRR22963925
rm -r SRR22963926
rm -r SRR22963927
rm -r SRR22963928
zcat SRR22963913_2.fastq.gz SRR22963914_2.fastq.gz SRR22963916_2.fastq.gz SRR22963917_2.fastq.gz SRR22963918_2.fastq.gz SRR22963919_2.fastq.gz SRR22963920_2.fastq.gz SRR22963921_2.fastq.gz SRR22963922_2.fastq.gz SRR22963923_2.fastq.gz SRR22963924_2.fastq.gz SRR22963925_2.fastq.gz SRR22963926_2.fastq.gz SRR22963927_2.fastq.gz SRR22963928_2.fastq.gz | gzip > SRS16347195.1.fq.gz
zcat SRR22963913_3.fastq.gz SRR22963914_3.fastq.gz SRR22963916_3.fastq.gz SRR22963917_3.fastq.gz SRR22963918_3.fastq.gz SRR22963919_3.fastq.gz SRR22963920_3.fastq.gz SRR22963921_3.fastq.gz SRR22963922_3.fastq.gz SRR22963923_3.fastq.gz SRR22963924_3.fastq.gz SRR22963925_3.fastq.gz SRR22963926_3.fastq.gz SRR22963927_3.fastq.gz SRR22963928_3.fastq.gz | gzip > SRS16347195.2.fq.gz
rm SRR22963913_1.fastq.gz SRR22963914_1.fastq.gz SRR22963916_1.fastq.gz SRR22963917_1.fastq.gz SRR22963918_1.fastq.gz SRR22963919_1.fastq.gz SRR22963920_1.fastq.gz SRR22963921_1.fastq.gz SRR22963922_1.fastq.gz SRR22963923_1.fastq.gz SRR22963924_1.fastq.gz SRR22963925_1.fastq.gz SRR22963926_1.fastq.gz SRR22963927_1.fastq.gz SRR22963928_1.fastq.gz
rm SRR22963913_2.fastq.gz SRR22963914_2.fastq.gz SRR22963916_2.fastq.gz SRR22963917_2.fastq.gz SRR22963918_2.fastq.gz SRR22963919_2.fastq.gz SRR22963920_2.fastq.gz SRR22963921_2.fastq.gz SRR22963922_2.fastq.gz SRR22963923_2.fastq.gz SRR22963924_2.fastq.gz SRR22963925_2.fastq.gz SRR22963926_2.fastq.gz SRR22963927_2.fastq.gz SRR22963928_2.fastq.gz
rm SRR22963913_3.fastq.gz SRR22963914_3.fastq.gz SRR22963916_3.fastq.gz SRR22963917_3.fastq.gz SRR22963918_3.fastq.gz SRR22963919_3.fastq.gz SRR22963920_3.fastq.gz SRR22963921_3.fastq.gz SRR22963922_3.fastq.gz SRR22963923_3.fastq.gz SRR22963924_3.fastq.gz SRR22963925_3.fastq.gz SRR22963926_3.fastq.gz SRR22963927_3.fastq.gz SRR22963928_3.fastq.gz

## for the BABOON data, we must download all samples manually

module add sratoolkit/3.0.0

# Baboon36506A

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Papio_anubis/SRS16347198 && cd $_
prefetch --max-size 100G SRR22963911
fastq-dump --split-files --gzip SRR22963911/SRR22963911.sra
rm -r SRR22963911
mv SRR22963911_2.fastq.gz SRS16347198.1.fq.gz
mv SRR22963911_3.fastq.gz SRS16347198.2.fq.gz
rm SRR22963911_1.fastq.gz

# Baboon36506B

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Papio_anubis/SRS16347196 && cd $_
prefetch --max-size 100G SRR22963915
fastq-dump --split-files --gzip SRR22963915/SRR22963915.sra
rm -r SRR22963915
mv SRR22963915_2.fastq.gz SRS16347196.1.fq.gz
mv SRR22963915_3.fastq.gz SRS16347196.2.fq.gz
rm SRR22963915_1.fastq.gz

# Baboon37979

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Papio_anubis/SRS16347197 && cd $_
prefetch --max-size 100G SRR22963912
fastq-dump --split-files --gzip SRR22963912/SRR22963912.sra
rm -r SRR22963912
mv SRR22963912_2.fastq.gz SRS16347197.1.fq.gz
mv SRR22963912_3.fastq.gz SRS16347197.2.fq.gz
rm SRR22963912_1.fastq.gz

## for the YAK data, we must download the sample manually

module add sratoolkit/3.0.0

mkdir /project/GorielyLab2021/sbush/ssc_atlas/fq/Bos_grunniens/SRS16827904 && cd $_
prefetch --max-size 100G SRR23546931 SRR23546932 SRR23546933 SRR23546934
fastq-dump --split-files --gzip SRR23546931/SRR23546931.sra
fastq-dump --split-files --gzip SRR23546932/SRR23546932.sra
fastq-dump --split-files --gzip SRR23546933/SRR23546933.sra
fastq-dump --split-files --gzip SRR23546934/SRR23546934.sra
rm -r SRR23546931
rm -r SRR23546932
rm -r SRR23546933
rm -r SRR23546934
zcat SRR23546931_1.fastq.gz SRR23546932_1.fastq.gz SRR23546933_1.fastq.gz SRR23546934_1.fastq.gz | gzip > SRS16827904.1.fq.gz
zcat SRR23546931_2.fastq.gz SRR23546932_2.fastq.gz SRR23546933_2.fastq.gz SRR23546934_2.fastq.gz | gzip > SRS16827904.2.fq.gz 
rm SRR23546931_1.fastq.gz SRR23546932_1.fastq.gz SRR23546933_1.fastq.gz SRR23546934_1.fastq.gz SRR23546931_2.fastq.gz SRR23546932_2.fastq.gz SRR23546933_2.fastq.gz SRR23546934_2.fastq.gz

=cut

use strict;
use warnings;

# REQUIREMENTS
my $root      = '/project/GorielyLab2021/sbush/ssc_atlas';
my $species   = 'Homo_sapiens'; # 'Mus_musculus'; # 'Macaca_mulatta'; # 'Gallus_gallus'; # 'Bubalus_bubalis'; # 'Rattus_norvegicus'; # 'Sus_scrofa'; # 'Ovis_aries'; # 'Bos_grunniens';
my $metadata  = "$root/prerequisites/metadata.$species.tsv"; # manually created
my $bam2fastq = '/project/GorielyLab2021/sbush/programs/bamtofastq-1.3.2';
my $ssh_file  = '/home/s/sbush/.aspera/connect/etc/asperaweb_id_dsa.openssh'; # see https://www.biostars.org/p/93482/
my $fatal     = 0;
if (!(-e($ssh_file)))  { $fatal++; print "ERROR: cannot find $ssh_file\n";  }
if (!(-e($metadata)))  { $fatal++; print "ERROR: cannot find $metadata\n";  }
if (!(-e($bam2fastq))) { $fatal++; print "ERROR: cannot find $bam2fastq\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_procs = 10;
my $ascp_or_wget = 'wget';

# OUTPUT
my $out_dir = "$root/fq";
if (!(-d($out_dir))) 			{ mkdir  $out_dir 			or die $!; }
if (!(-d("$out_dir/$species"))) { mkdir "$out_dir/$species" or die $!; }
my $out_sh  = "$root/run_fq_downloader.$species.sh";
open(SH,'>',$out_sh) or die $!;
print SH "#!/bin/bash\n";
print SH "module add aspera/3.9.8\n";

# WHAT SAMPLE IDs ARE WE GOING TO RUN?
my %samples = (); my %library_layouts = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $study = $line[0]; my $sample_id = $line[1]; my $run_ids = $line[2]; my $library_layout = $line[3]; my $techno = $line[6]; my $bam_urls = $line[7]; my $phenotype = $line[10];
	  my @run_ids = split(/\,/,$run_ids); my @bam_urls = split(/\,/,$bam_urls);
	  $samples{$sample_id}{layout} = $library_layout;
	  $samples{$sample_id}{techno} = $techno;
	  $library_layouts{$sample_id}{$library_layout}++;
	  if ($bam_urls ne 'NA') # [BAM ROUTE] ONLY APPLICABLE FOR 10X DATA THAT NEEDS TO BE UNPACKED
		{ for(my $x=0;$x<@run_ids;$x++)
			{ my $bam_url = $bam_urls[$x]; my $run_id = $run_ids[$x];
			  push(@{$samples{$sample_id}{bam_urls}},[$bam_url,$run_id]);
			}
		}
	  elsif ($bam_urls eq 'NA') # [FQ ROUTE]
		{ foreach my $run_id (@run_ids)
			{ my $first_3; my $first_6; my $digits;
			  if ($run_id =~ /^(.{3}).*?$/) { $first_3 = $1; }
			  if ($run_id =~ /^(.{6}).*?$/) { $first_6 = $1; }
			  if ($run_id =~ /^.*?(\d+)$/)  { $digits  = $1; }
			  my $number_of_digits = length($digits);
			  my $ena_url_1 = ''; my $ena_url_2 = '';
			  my $wget_url_1 = ''; my $wget_url_2 = '';
			  if ($number_of_digits == 6)
				{ if ($library_layout eq 'PAIRED')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
					  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id"."_1.fastq.gz";
					  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id"."_2.fastq.gz";
					}
				  elsif ($library_layout eq 'SINGLE')
					{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$run_id/$run_id.fastq.gz";
					  $ena_url_2  = 'NA';
					  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$run_id/$run_id.fastq.gz";
					  $wget_url_2 = 'NA';
					}
				}
			  elsif ($number_of_digits == 7)
				{ if ($digits =~ /^.+?(\d{1})$/)
					{ my $last_digit = $1;
					  if ($library_layout eq 'PAIRED')
						{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
						  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
						  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_1.fastq.gz";
						  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id"."_2.fastq.gz";
						}
					  elsif ($library_layout eq 'SINGLE')
						{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id.fastq.gz";
						  $ena_url_2  = 'NA';
						  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/00$last_digit/$run_id/$run_id.fastq.gz";
						  $wget_url_2 = 'NA';
						}
					}
				}
			  elsif ($number_of_digits == 8)
				{ if ($digits =~ /^.+?(\d{2})$/)
					{ my $last_two_digits = $1;
					  if ($library_layout eq 'PAIRED')
						{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
						  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
						  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_1.fastq.gz";
						  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id"."_2.fastq.gz";
						}
					  elsif ($library_layout eq 'SINGLE')
						{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id.fastq.gz";
						  $ena_url_2  = 'NA';
						  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/0$last_two_digits/$run_id/$run_id.fastq.gz";
						  $wget_url_2 = 'NA';
						}
					}
				}
			  elsif ($number_of_digits == 9)
				{ if ($digits =~ /^.+?(\d{3})$/)
					{ my $last_three_digits = $1;
					  if ($library_layout eq 'PAIRED')
						{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
						  $ena_url_2  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
						  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_1.fastq.gz";
						  $wget_url_2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id"."_2.fastq.gz";
						}
					  elsif ($library_layout eq 'SINGLE')
						{ $ena_url_1  = "era-fasp\@fasp.sra.ebi.ac.uk:/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id.fastq.gz";
						  $ena_url_2  = 'NA';
						  $wget_url_1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$first_6/$last_three_digits/$run_id/$run_id.fastq.gz";
						  $wget_url_2 = 'NA';
						}
					}
				}
			  push(@{$samples{$sample_id}{fq_urls}},[$ena_url_1,$ena_url_2,$wget_url_1,$wget_url_2]);
			}
		}
	}
close(IN) or die $!;

# FOR EACH SAMPLE...
my @samples = ();
while((my $sample_id,my $irrel)=each(%samples))
	{ push(@samples,$sample_id); }
my @sorted_samples = sort {$a cmp $b} @samples;
my $files_seen = 0; my $files_total = @sorted_samples;
foreach my $sample_id (@sorted_samples)
	{ $files_seen++;
	  print "$files_seen of $files_total\n";
	  my $num_library_layouts = scalar keys %{$library_layouts{$sample_id}};
	  next if ($num_library_layouts != 1); # CHECKPOINT: skip if the sample ID is associated with both paired-end and single-end fqs; these have already been pre-processed in some way
	  my $library_layout = '';
	  while((my $this_layout,my $irrel)=each(%{$library_layouts{$sample_id}}))
		{ $library_layout = $this_layout; }
	  
	  # SKIP IF WE'VE OBTAINED THESE FQs ALREADY
	  next if ( ($library_layout eq 'SINGLE') and (-e("$out_dir/$species/$sample_id/$sample_id.fq.gz")) );
	  next if ( ($library_layout eq 'PAIRED') and (-e("$out_dir/$species/$sample_id/$sample_id.1.fq.gz")) and (-e("$out_dir/$species/$sample_id/$sample_id.2.fq.gz")) );

	  # IF WE HAVEN'T OBTAINED THESE FQs ALREADY BUT WE HAVE CREATED THE OUTPUT DIR (PERHAPS IN A PREVIOUS FAILED RUN?), THEN PURGE ITS CONTENTS - WE DON'T TRUST THEM
	  if (-e("$out_dir/$species/$sample_id"))
		{ print SH "rm -r $out_dir/$species/$sample_id\n"; }
	  print SH "mkdir $out_dir/$species/$sample_id\n";
	  print SH "cd $out_dir/$species/$sample_id\n";
	  
	  my $fq1 = "$out_dir/$species/$sample_id/$sample_id.1.fq.gz";
	  my $fq2 = "$out_dir/$species/$sample_id/$sample_id.2.fq.gz";
	  my $fq  = "$out_dir/$species/$sample_id/$sample_id.fq.gz";
	  if (exists($samples{$sample_id}{fq_urls}))
		{ # [FQ ROUTE] DOWNLOAD ALL FQs FOR EACH RUN OF THIS SAMPLE. THE VAST MAJORITY OF SAMPLES IN $metadata WILL TAKE THIS ROUTE.
		  my $file_line1 = ''; my $file_line2 = '';
		  my @fqs_to_delete = ();
		  my @arr = @{$samples{$sample_id}{fq_urls}};
		  for(my $x=0;$x<@arr;$x++)
			{ my $url1 = $arr[$x][0]; my $url2 = $arr[$x][1]; my $wget1 = $arr[$x][2]; my $wget2 = $arr[$x][3];
			  my $fq_1 = $url1; my $fq_2 = $url2;
			  if ($fq_1 =~ /^.+\/(.*?)$/) { $fq_1 = $1; }
			  if ($fq_2 =~ /^.+\/(.*?)$/) { $fq_2 = $1; }
			  if ($ascp_or_wget eq 'ascp')
				{ print SH "ascp -QT -l 100m -P33001 -i $ssh_file $url1 $fq_1\n"; }
			  elsif ($ascp_or_wget eq 'wget')
				{ print SH "wget $wget1 -O $fq_1\n"; }
			  if (($url2 eq 'NA') and ($library_layout eq 'SINGLE'))
				{ $file_line1 .= "$out_dir/$species/$sample_id/$fq_1 ";
				  push(@fqs_to_delete,$fq_1);
				}
			  elsif (($url2 ne 'NA') and ($library_layout eq 'PAIRED'))
				{ if ($ascp_or_wget eq 'ascp')
					{ print SH "ascp -QT -l 100m -P33001 -i $ssh_file $url2 $fq_2\n"; }
				  elsif ($ascp_or_wget eq 'wget')
					{ print SH "wget $wget2 -O $fq_2\n"; }
				  $file_line1 .= "$out_dir/$species/$sample_id/$fq_1 ";
				  $file_line2 .= "$out_dir/$species/$sample_id/$fq_2 ";
				  push(@fqs_to_delete,$fq_1,$fq_2);
				}
			}
		
		  # [FQ ROUTE] IF THERE HAVE BEEN MULTIPLE SETS OF RUN FQS PER SAMPLE, THEN MERGE THEM ALL TOGETHER INTO ONE FINAL SAMPLE.FQ
		  if ($#arr > 0)
			{ if ($library_layout eq 'PAIRED')
				{ print SH "cat $file_line1 > $fq1\n";
				  print SH "cat $file_line2 > $fq2\n";
				}
			  elsif ($library_layout eq 'SINGLE')
				{ print SH "cat $file_line1 > $fq\n";
				}
			  foreach my $fq (@fqs_to_delete)
				{ print SH "rm $fq\n"; }
			}
		  elsif ($#arr == 0) # ELSE JUST RENAME THE SINGLE RUN.FQ TO SAMPLE.FQ
			{ if ($library_layout eq 'PAIRED')
				{ print SH "mv $file_line1 $fq1\n";
				  print SH "mv $file_line2 $fq2\n";
				}
			  elsif ($library_layout eq 'SINGLE')
				{ print SH "mv $file_line1 $fq\n";
				}
			}
		}
	  elsif (exists($samples{$sample_id}{bam_urls}))
		{ # [BAM ROUTE] DOWNLOAD THE 10X BAMS DIRECTLY AND THEN UNPACK TO FQ
		  my @arr = @{$samples{$sample_id}{bam_urls}};
		  for(my $x=0;$x<@arr;$x++)
			{ my $url = $arr[$x][0]; my $run_id = $arr[$x][1];
			  print SH "wget -c \"$url\" -O $run_id.bam\n";
			  print SH "$bam2fastq --nthreads $num_procs --reads-per-fastq 10000000000000000 $run_id.bam $run_id-fqs\n";
			  print SH "cat $run_id-fqs/*/*_R1_001.fastq.gz > $fq1\n";
			  print SH "cat $run_id-fqs/*/*_R2_001.fastq.gz > $fq2\n";
			  print SH "rm -r $run_id-fqs\n";
			}
		}
	}
close(SH) or die $!;
exit 1;