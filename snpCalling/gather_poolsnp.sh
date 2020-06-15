#!/usr/bin/env bash
#
#SBATCH -J split_and_run # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch ${wd}/DEST/snpCalling/gather_poolsnp.sh
### sacct -j 12825841
### cat /scratch/aob2x/dest/slurmOutput/split_and_run.12825614
module load htslib bcftools intel/18.0 intelmpi/18.0 R/3.6.0


wd="/scratch/aob2x/dest"
outdir="/scratch/aob2x/dest/sub_vcfs"

#head -n5 ${wd}/dest/poolSNP_jobs.csv | sed 's/,/_/g' | sed 's/^/\/scratch\/aob2x\/dest\/sub_vcfs\//g' | sed 's/$/.bcf/g' > /scratch/aob2x/dest/sub_vcfs/bcfs_order
#ls -d $outdir/*.vcf.gz > /scratch/aob2x/dest/sub_vcfs/vcfs_order
# cat /scratch/aob2x/dest/sub_vcfs/vcfs_order | sort -t"_" -k2,2 -k3g,3  > /scratch/aob2x/dest/sub_vcfs/vcfs_order.sort
#
 bcftools concat \
 --threads 20 \
 -f /scratch/aob2x/dest/sub_vcfs/vcfs_order.sort \
 -O v \
 -n \
 -o /scratch/aob2x/dest/dest.June14_2020.bcf

#mv /scratch/aob2x/dest/dest.June14_2020.vcf /scratch/aob2x/dest/dest.June14_2020.bcf

bcftools view /scratch/aob2x/dest/dest.June14_2020.bcf > /scratch/aob2x/dest/dest.June14_2020.vcf

java -jar ~/snpEff/snpEff.jar eff BDGP6.86 /scratch/aob2x/dest/dest.June14_2020.vcf -v > /scratch/aob2x/dest/dest.June14_2020.ann.vcf

Rscript --vanilla ${wd}/DEST/snpCalling/vcf2gds.R
