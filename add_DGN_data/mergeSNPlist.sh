#!/usr/bin/env bash
#
#SBATCH -J mergeSNPlist # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/mergeSNPlist.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/mergeSNPlist.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/7.1.0 bedops/2.4.1

bedops -u \
/scratch/aob2x/dest/dgn/sitesData/dgn_sites.dm6.bed \
/scratch/aob2x/dest/drosRTEC/drosRTEC_sites.dm6.bed \
/scratch/aob2x/dest/drosEU/drosEU_sites.dm6.bed |
grep -E "chr2L|chr2R|chr3L|chr3R|chrX" > /scratch/aob2x/dest/dest/dgn_drosRTEC_drosEU.sites.dm6.bed
