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

### 2a. Make PoolSNP based VCF file (bgzip out)
```bash
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh 001
sbatch --array=1-$( wc -l ${wd}/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/snpCalling/run_poolsnp.sh 01
```

### 2b. Collect PoolSNP (bcf out)
```bash
sbatch --array=1-7 ${wd}/DEST/snpCalling/gather_poolsnp.sh 001
sbatch --array=1-7 ${wd}/DEST/snpCalling/gather_poolsnp.sh 01
```

### 3c. Bind chromosomes, annotate and convert (bgzip out; GDS out)
```bash
sbatch ${wd}/DEST/snpCalling/annotate.sh 001
sbatch ${wd}/DEST/snpCalling/annotate.sh 01

```









## SNP calling with PoolSNP based on SYNC file format

call SNPs on joined SYNC file based on (1) minimum allele count (across all samples) or (2) minimum allele frequency (across all samples). Note, that only positions with a proportion of missing data <= than the miss-frac threshold are retained. The max-cov argument is only used for the VCF header.

```bash
python /Users/mkapun/Documents/GitHub/DEST/PoolSNP4Sync/PoolSnp.py \
--sync testMpileup2Sync/joined_masked.sync.gz \
--min-cov 4 \
--max-cov 0.7 \
--min-count 4 \
--min-freq 0.01 \
--miss-frac 0.5 \
--names sample1,sample2 \
> testMpileup2Sync/joined.vcf
```

## make job file to split out chromosomal sub chunks
  ### makes job DEST/PoolSNP4Sync/makeJobs.R

  ### produces /scratch/aob2x/dest/dest/poolSNP_jobs.csv
