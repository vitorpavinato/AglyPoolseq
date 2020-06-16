#!/usr/bin/env bash
#
#SBATCH -J annotate # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

echo ${1}

#bcftools view /scratch/aob2x/dest/dest.June14_2020.bcf > /scratch/aob2x/dest/dest.June14_2020.vcf
#
#java -jar ~/snpEff/snpEff.jar eff BDGP6.86 /scratch/aob2x/dest/dest.June14_2020.vcf -v > /scratch/aob2x/dest/dest.June14_2020.ann.vcf
#
#Rscript --vanilla ${wd}/DEST/snpCalling/vcf2gds.R
