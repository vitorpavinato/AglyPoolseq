#!/usr/bin/env bash
#
#SBATCH -J download_DGN # A single job name for the array
#SBATCH --ntasks-per-node=1 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 48 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/download_DGN.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/download_DGN.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
