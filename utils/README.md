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
  #sacct -j 16534019 | grep "TIMEOUT" | cut -f1 -d' ' | cut -f2 -d'_' | tr '\n' ','


  sbatch --array=2,4,5,6,7,8,9,11,14,15,16,17,18,19,20,21,22,23,25,26,28,29,30,31,32,33,34,37,38,39,40,41,43,46,47,49,50,51,52,55,56,57,58,59,60,61,62,63,64,65,66,69,70,71,72,74,75,76,77,78,79,81,82,84,85,87,88,89,90,91,93,94,95,98,99,100,102,104,110,113,114,115,117,118,119,120,121,122,123,124,126,128,139,142,143,144,145,150,151,152,153,154,155,157,161,163,170,173,174,177,182,183,188,192,195,196,197,200,203,206,210,211,215,217,218,224,236,237,240,242,246 \
  ${wd}/DEST/utils/getPrivateSNPs.sh
  sacct -j 16901160

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
