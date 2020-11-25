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
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh all PoolSNP 001 50 10Nov2020 poolSNP_jobs.csv
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh PoolSeq PoolSNP 001 50 10Nov2020 poolSNP_jobs.csv
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh PoolSeq SNAPE NA NA 10Nov2020 poolSNP_jobs.csv
```

### 2b. Collect PoolSNP (bcf out)
```bash
sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh all PoolSNP 001 50 10Nov2020
sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh PoolSeq PoolSNP 001 50 10Nov2020
sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh PoolSeq SNAPE NA NA 10Nov2020
```


### 2c. Bind chromosomes, annotate and convert (bgzip out; GDS out)
```bash
sbatch ${wd}/DEST/snpCalling/annotate.sh all PoolSNP 001 50 10Nov2020
sbatch ${wd}/DEST/snpCalling/annotate.sh PoolSeq PoolSNP 001 50 10Nov2020
sbatch ${wd}/DEST/snpCalling/annotate.sh PoolSeq SNAPE NA NA 10Nov2020
```



## 3. Parameter evaluation for PoolSNP (global MAC & MAF thresholds)
### 3a. Random sample of ~10% of data:
```bash
  shuf -n 100 ${wd}/poolSNP_jobs.csv > ${wd}/poolSNP_jobs.sample.csv
```

### 3b. Run pool_snp
```bash
  module load parallel

  runJob () {
    wd="/scratch/aob2x/dest"
    sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.sample.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh all PoolSNP ${1} ${2} paramTest poolSNP_jobs.sample.csv
  }
  export -f runJob

  parallel -j 1 runJob ::: 001 01 05 ::: 5 10 15 20 50 100
```
sacct -j 19030434
sacct -j 19030435
sacct -j 19030436
sacct -j 19030437
sacct -j 19030438
sacct -j 19030439
sacct -j 19030440
sacct -j 19030441
sacct -j 19030442
sacct -j 19030443
sacct -j 19030444
sacct -j 19030445
sacct -j 19030446
sacct -j 19030447
sacct -j 19030448
sacct -j 19030449
sacct -j 19030450
sacct -j 19030451


### 3c. Run modified gather
```bash
  module load parallel

  runJob () {
    wd="/scratch/aob2x/dest"
    sbatch --array=1-8 ${wd}/DEST/snpCalling/gather_poolsnp.sh all PoolSNP ${1} ${2} paramTest
  }
  export -f runJob

  parallel -j 1 runJob ::: 001 01 05 ::: 5 10 15 20 50 100
```
sacct -j 19029899
sacct -j 19029900
sacct -j 19029901
sacct -j 19029902
sacct -j 19029903
sacct -j 19029904
sacct -j 19029905
sacct -j 19029906
sacct -j 19029907
sacct -j 19029908
sacct -j 19029909
sacct -j 19029910
sacct -j 19029911
sacct -j 19029912
sacct -j 19029913
sacct -j 19029914
sacct -j 19029915
sacct -j 19029916






### 3d. Bind, annotate and convert (bgzip out; GDS out)
```bash
sbatch --array=1-$( ls /scratch/aob2x/dest/sub_bcf_paramTest/dest.June14_2020.maf001.*.paramTest.bcf | rev | cut -d'.' -f3,4 | rev | sort | uniq | awk '{print NR}' | tail -n1 ) ${wd}/DEST/snpCalling/annotate_paramtest.sh
```


## 4. Generate liftover table for polymorphic sites (dm6 -> dm3/r5)
