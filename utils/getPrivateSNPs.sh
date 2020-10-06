#!/usr/bin/env bash
#
#SBATCH -J getPrivate # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/getPrivate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/getPrivate.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab



### run as: sbatch --array=1-246 ${wd}/DEST/utils/getPrivateSNPs.sh
### sacct -j 13029741

module load htslib bcftools intel/18.0 intelmpi/18.0 parallel

wd="/scratch/aob2x/dest"

#SLURM_ARRAY_TASK_ID=1


Rscript ${wd}/DEST/utils/getPrivateSNPs.R ${SLURM_ARRAY_TASK_ID} PoolSNP
Rscript ${wd}/DEST/utils/getPrivateSNPs.R ${SLURM_ARRAY_TASK_ID} SNAPE
