#!/usr/bin/env bash
#
#SBATCH -J getSNP_list_drosRTEC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/getSNP_list_drosRTEC.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/getSNP_list_drosRTEC.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#SLURM_ARRAY_TASK_ID=2



  ### in bash
/mnt/spicy_1/sin/liftOver \
/mnt/spicy_2/dest/drosRTEC_sites.dm3.bed \
/mnt/spicy_1/sin/dm3ToDm6.over.chain \
/mnt/spicy_2/dest/drosRTEC_sites.dm6.bed \
/mnt/spicy_2/dest/drosRTEC_sites.unmapped.bed
