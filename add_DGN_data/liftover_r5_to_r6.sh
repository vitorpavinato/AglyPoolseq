#!/usr/bin/env bash
#
#SBATCH -J liftover_r5_to_r6 # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:30:00 ### 0.5 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/liftover_r5_to_r6.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/liftover_r5_to_r6.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#SLURM_ARRAY_TASK_ID=2

if [ ${SLURM_ARRAY_TASK_ID} -eq 1 ]; then
  chr=2L
elif [ ${SLURM_ARRAY_TASK_ID} -eq 2 ]; then
  chr=2R
elif [ ${SLURM_ARRAY_TASK_ID} -eq 3 ]; then
  chr=3L
elif [ ${SLURM_ARRAY_TASK_ID} -eq 4 ]; then
  chr=3R
elif [ ${SLURM_ARRAY_TASK_ID} -eq 5 ]; then
  chr=X
fi

echo $chr

~/liftOver \
/scratch/aob2x/dest/dgn/sitesData/dgn_sites_${chr}.bed \
/scratch/aob2x/dest/dgn/liftoverChains/dm3ToDm6.over.chain \
/scratch/aob2x/dest/dgn/sitesData/dgn_sites_${chr}.dm6.bed \
/scratch/aob2x/dest/dgn/sitesData/dgn_sites_${chr}.unmapped.bed
