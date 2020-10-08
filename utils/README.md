# Some scripts to generate shared data files after SNP calling

## 0. Define working directory
```bash
wd=/scratch/aob2x/dest
```

## 1. bamMappingRates.R
  ### makes mps? need to check

## 2. Get counts of private polymorphisms using only PoolSeq samples for both PoolSNP & SNAPE.
### 2a. Generate output for each population
```bash
sbatch --array=1-246 ${wd}/DEST/utils/getPrivateSNPs.sh
```
#sbatch --array=1 ${wd}/DEST/utils/getPrivateSNPs.sh
#sacct -j 16534019
