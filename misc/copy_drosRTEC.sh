#!/usr/bin/env bash
#
#SBATCH -J copy_drosRTEC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 48:00:00 ### 48 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/copy_drosRTEC.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/copy_drosRTEC.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

cp -avr /nv/vol186/bergland-lab/alan/mapped /scratch/aob2x/dest/drosRTEC/raw_n_mapped_data/.
