#!/usr/bin/env bash
#
#SBATCH -J cp_nv2scratch # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 165:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/cp_nv2scratch.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/cp_nv2scratch.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

cp -vr /nv/vol186/bergland-lab/alan/mapped /project/berglandlab/alan/drosRTEC/.
