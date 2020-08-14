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

### 2a. Make PoolSNP based VCF file (bgzip out). Uses MAF > 0.001 & MAC > 50. These are reasonable thresholds that produce consistent pn/ps, number of SNPs, et, but can be filtered at a later stage using standard VCF tools.
```bash
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh 001 50
```
sacct -j 14256378

### 2b. Collect PoolSNP (bcf out)
```bash
sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh 001 50
```


### 2c. Bind chromosomes, annotate and convert (bgzip out; GDS out)
```bash
sbatch ${wd}/DEST/snpCalling/annotate.sh 001 50 Aug9_2020
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
