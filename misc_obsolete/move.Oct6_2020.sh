#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 0.5 hours
##SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.16425773.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

### sbatch ${wd}/DEST/misc_obsolete/move.Oct6_2020.sh
# sacct -j 16428909

wd="/scratch/aob2x/dest"
rsync ${wd}/dest.PoolSeq.SNAPE.001.50.ann.gds /project/berglandlab/DEST/gds/.
rsync ${wd}/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz /project/berglandlab/DEST/vcf/.
rsync ${wd}/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz.tbi /project/berglandlab/DEST/vcf/.
rsync ${wd}/dest.PoolSeq.SNAPE.001.50.header.bcf /project/berglandlab/DEST/bcf/.
