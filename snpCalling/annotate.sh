#!/usr/bin/env bash
#
#SBATCH -J annotate # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 6 hours
#SBATCH --mem 90G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

maf=${1}
#maf=001
wd=/scratch/aob2x/dest


bcftools concat \
${wd}/sub_bcf/dest.June14_2020.maf001.2L.${maf}.bcf \
${wd}/sub_bcf/dest.June14_2020.maf001.2R.${maf}.bcf \
${wd}/sub_bcf/dest.June14_2020.maf001.3L.${maf}.bcf \
${wd}/sub_bcf/dest.June14_2020.maf001.3R.${maf}.bcf \
${wd}/sub_bcf/dest.June14_2020.maf001.X.${maf}.bcf \
${wd}/sub_bcf/dest.June14_2020.maf001.4.${maf}.bcf \
${wd}/sub_bcf/dest.June14_2020.maf001.Y.${maf}.bcf \
-n \
-o ${wd}/dest.June14_2020.${maf}.bcf

bcftools view \
--threads 10 \
${wd}/dest.June14_2020.${maf}.bcf > ${wd}/dest.June14_2020.${maf}.vcf

java -jar ~/snpEff/snpEff.jar \
eff \
-v \
-t 10 \
BDGP6.86 \
${wd}/dest.June14_2020.${maf}.vcf > \
${wd}/dest.June14_2020.${maf}.ann.vcf

Rscript --vanilla ${wd}/DEST/snpCalling/vcf2gds.R ${wd}/dest.June14_2020.${maf}.ann.vcf

bgzip -c ${wd}/dest.June14_2020.${maf}.ann.vcf > ${wd}/dest.June14_2020.${maf}.ann.vcf.gz
