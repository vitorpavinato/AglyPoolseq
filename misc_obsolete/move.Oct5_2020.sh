#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 0.5 hours
##SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

### sbatch ${wd}/DEST/misc_obsolete/move.Aug22_2020.001.50.sh
# sacct -j 14487964

wd="/scratch/aob2x/dest"
rsync ${wd}/dest.PoolSeq.PoolSNP.001.50.ann.gds /project/berglandlab/DEST/
rsync ${wd}/dest.PoolSeq.PoolSNP.001.50.ann.vcf.gz /project/berglandlab/DEST/
rsync ${wd}/dest.PoolSeq.PoolSNP.001.50.ann.vcf.gz.tbi /project/berglandlab/DEST/
#rsync ${wd}/dest.Aug9_2020.001.50.header.bcf /project/berglandlab/DEST/dest.July6_2020.001.10.bcf
