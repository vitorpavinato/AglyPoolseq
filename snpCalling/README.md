# Scripts to call polymorphic sites using PoolSNP and SNAPE.

## Description
> This set of scripts generates VCF files from whole-genome SYNC files using one of two SNP calling pipelines.

### 0. Define working directory
```bash
wd=/scratch/aob2x/dest
```

### 1. Generate job id file
```bash
Rscript ${wd}/DEST/snpCalling/makeJobs.R
```

### 2a. Make PoolSNP based VCF file (bgzip out). Uses MAF > 0.001 & MAC > 50. These are reasonable thresholds that produce consistent pn/ps, number of SNPs, et, but can be filtered at a later stage using standard VCF tools. </br>
First paramter is the population set ('all' samples or just the 'PoolSeq' samples). Second parameter is the SNP calling method (PoolSNP or SNAPE). If method == PoolSNP, third parameter is MAF filter, fourth is MAC filter. These are retained for the SNAPE version just to keep things consistent.

```bash
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh all PoolSNP 001 50
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh PoolSeq PoolSNP 001 50
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh PoolSeq SNAPE 001 50
```
sbatch --array=1 ${wd}/DEST/snpCalling/run_poolsnp.sh PoolSeq PoolSNP 001 50
sbatch --array=1 ${wd}/DEST/snpCalling/run_poolsnp.sh PoolSeq SNAPE 001 50

sacct -j 16292561
sacct -j 16292723

ls -lh 2L_1_138315*

cat /scratch/aob2x/dest/slurmOutput/split_and_run/

### 2b. Collect PoolSNP (bcf out)
```bash
#sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh 001 50 Aug22_2020
sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh PoolSeq PoolSNP 001 50
sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh PoolSeq SNAPE 001 50
```

sacct -j 16274411
sacct -j 16274414

### 2c. Bind chromosomes, annotate and convert (bgzip out; GDS out)
```bash
sbatch ${wd}/DEST/snpCalling/annotate.sh 001 50 Aug22_2020
```



## 3. Parameter evaluation. Some modifications to introduce internal parallel function. Runs through 3c, I think.
### 3a. Random sample of ~10% of data:
```bash
shuf -n 100 ${wd}/poolSNP_jobs.csv > ${wd}/poolSNP_jobs.sample.csv
```
### 3b. Run modified pool_snp
```bash
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.sample.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp_paramtest.sh
```
### 3c. Run modified gather
```bash
sbatch --array=1-$( cat ${wd}/poolSNP_jobs.sample.csv | cut -f1 -d',' | sort | uniq | awk '{print NR}' | tail -n1 ) ${wd}/DEST/snpCalling/ gather_poolsnp_paramtest.sh
```
### 3d. Bind, annotate and convert (bgzip out; GDS out)
```bash
sbatch --array=1-$( ls /scratch/aob2x/dest/sub_bcf_paramTest/dest.June14_2020.maf001.*.paramTest.bcf | rev | cut -d'.' -f3,4 | rev | sort | uniq | awk '{print NR}' | tail -n1 ) ${wd}/DEST/snpCalling/annotate_paramtest.sh
```


## 4.
