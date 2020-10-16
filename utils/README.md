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

  tar czvf /project/berglandlab/DEST/privateSNPs/privateSNPs.SNAPE.tar.gz ${wd}/privateSNPs/*SNAPE.csv
  tar czvf /project/berglandlab/DEST/privateSNPs/privateSNPs.PoolSNP.tar.gz ${wd}/privateSNPs/*PoolSNP.csv

  tar czvf /project/berglandlab/DEST/privateSNPs/privateSNPs.summary.tar.gz ${wd}/privateSNPs/*summary.delim
  
  ```
  sacct -j 16914500

## 3. Make liftover of polymorphic locations to dm3 coordinates
### 3a. download liftover file; liftOver utility assumed to exist somewhere
  ```bash
  wget \
  -O ${wd}/referenceGenome/dm6ToDm3.over.chain.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/dm6/liftOver/dm6ToDm3.over.chain.gz
  ```
### 3b. make bed file
  ```bash
  zcat ${wd}/dest.all.PoolSNP.001.50.ann.vcf.gz | grep -v "#" | cut -f1,2 |  awk '{print "chr"$1"\t"$2"\t"$2+1"\tdm6_"$1"_"$2}' > ${wd}/dest.all.PoolSNP.001.50.dm6.bed
  ```

### 3c. liftover to dm3 and parse
  ```bash
  ~/liftOver \
  ${wd}/dest.all.PoolSNP.001.50.dm6.bed \
  ${wd}/referenceGenome/dm6ToDm3.over.chain.gz \
  ${wd}/dest.all.PoolSNP.001.50.dm3.bed \
  ${wd}/dest.all.PoolSNP.001.50.dm6.unmapped

  sed 's/_/\t/g' ${wd}/dest.all.PoolSNP.001.50.dm3.bed | sed 's/chr//g' | awk '
  BEGIN {
    print "dm3_chr,dm3_pos,dm6_chr,dm6_pos"
  }
  {
    print $1","$2","$5","$6
  }
  ' | gzip -c - > ${wd}/DEST/utils/dest.all.PoolSNP.001.50.dm3.dm6.csv.gz
