#!/usr/bin/env bash
#
#SBATCH -J split_and_run # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch --array=1-$( wc -l ${wd}/dest/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/PoolSNP4Sync/run_poolsnp.sh
### sbatch --array=1-5 ${wd}/DEST/PoolSNP4Sync/run_poolsnp.sh
### sacct -j 12755543
module load htslib bcftools

wd="/scratch/aob2x/dest"
outdir="/scratch/aob2x/dest/sub_vcfs"

#head -n5 ${wd}/dest/poolSNP_jobs.csv | sed 's/,/_/g' | sed 's/^/\/scratch\/aob2x\/dest\/sub_vcfs\//g' | sed 's/$/.bcf/g' > /scratch/aob2x/dest/sub_vcfs/bcfs_order
ls -d $outdir/*.bcf > /scratch/aob2x/dest/sub_vcfs/bcfs_order

bcftools concat \
-f /scratch/aob2x/dest/sub_vcfs/bcfs_order \
-O v \
-o ~/test.all.vcf

bcftools view ~/test.all.vcf | less -S

grep "5735" ~/test.all.vcf | tr '\t' '\n'
