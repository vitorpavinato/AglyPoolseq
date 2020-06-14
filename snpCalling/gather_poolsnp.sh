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

### run as: sbatch -${wd}/DEST/PoolSNP4Sync/gather_poolsnp.sh
### sacct -j 12755543
module load htslib bcftools intel/18.0 intelmpi/18.0 R/3.6.0


wd="/scratch/aob2x/dest"
outdir="/scratch/aob2x/dest/sub_vcfs"

#head -n5 ${wd}/dest/poolSNP_jobs.csv | sed 's/,/_/g' | sed 's/^/\/scratch\/aob2x\/dest\/sub_vcfs\//g' | sed 's/$/.bcf/g' > /scratch/aob2x/dest/sub_vcfs/bcfs_order
ls -d $outdir/*.vcf.gz > /scratch/aob2x/dest/sub_vcfs/vcfs_order

bcftools concat \
--threads 20 \
-a \
-f /scratch/aob2x/dest/sub_vcfs/vcfs_order \
-O v \
-o /scratch/aob2x/dest.June14_2020.vcf

java -jar ~/snpEff/snpEff.jar eff BDGP6.86 /scratch/aob2x/dest.June14_2020.vcf > /scratch/aob2x/dest.June14_2020.ann.vcf

Rscript --vanilla ${wd}/DEST/snpCalling/vcf2gds.R
