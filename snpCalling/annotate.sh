#!/usr/bin/env bash
#
#SBATCH -J annotate # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 6 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load htslib bcftools intel/18.0 intelmpi/18.0 parallel R/3.6.3

maf=${1}
mac=${2}
stem=${3}

#maf=001;mac=50;stem=Aug9_2020

wd=/scratch/aob2x/dest

echo "concat"
  bcftools concat \
  ${wd}/sub_bcf/dest.${stem}.2L.${maf}.${mac}.bcf \
  ${wd}/sub_bcf/dest.${stem}.2R.${maf}.${mac}.bcf \
  ${wd}/sub_bcf/dest.${stem}.3L.${maf}.${mac}.bcf \
  ${wd}/sub_bcf/dest.${stem}.3R.${maf}.${mac}.bcf \
  ${wd}/sub_bcf/dest.${stem}.X.${maf}.${mac}.bcf \
  ${wd}/sub_bcf/dest.${stem}.4.${maf}.${mac}.bcf \
  ${wd}/sub_bcf/dest.${stem}.Y.${maf}.${mac}.bcf \
  -o ${wd}/dest.${stem}.${maf}.${mac}.bcf

echo "convert to vcf & annotate"
  bcftools view \
  --threads 10 \
  ${wd}/dest.${stem}.${maf}.${mac}.bcf | \
  java -jar ~/snpEff/snpEff.jar \
  eff \
  BDGP6.86 - > \
  ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf

echo "fix header" #this is now fixed in PoolSNP.py
  sed -i '0,/CHROM/{s/AF,Number=1/AF,Number=A/}' ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf
  sed -i '0,/CHROM/{s/AC,Number=1/AC,Number=A/}' ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf
  sed -i '0,/CHROM/{s/AD,Number=1/AD,Number=A/}' ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf
  sed -i '0,/CHROM/{s/FREQ,Number=1/FREQ,Number=A/}' ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf

  bcftools view -h ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf > ${wd}/tmp.header

  bcftools reheader --threads 10 -h ${wd}/tmp.header -o ${wd}/dest.${stem}.${maf}.${mac}.header.bcf ${wd}/dest.${stem}.${maf}.${mac}.bcf

echo "make GDS"
  Rscript --vanilla ${wd}/DEST/snpCalling/vcf2gds.R ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf

echo "bgzip & tabix"
  bgzip -c ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf > ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf.gz
  tabix -p vcf ${wd}/dest.${stem}.${maf}.${mac}.ann.vcf.gz
