#!/usr/bin/env bash
#
#SBATCH -J pnpsPrivate # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00 ### 6 hours
#SBATCH --mem 64G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/pnpsPrivate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/pnpsPrivate.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab



### run as: sbatch --array=2-246 ${wd}/DEST/utils/pnps.PrivateSNPs.sh
### sacct -j 17166145

### run as: sbatch --array=1 ${wd}/DEST/utils/pnps.PrivateSNPs.sh
### sacct -j 17164717

#sacct -j 17166145 | grep "OUT_OF" | grep -v "batch" | cut -f1 -d' ' | cut -f2 -d'_' | tr '\n' ','

### sbatch --array=3,4,7,28,75,80,81,102,114,141,151,152,157,167,168,185,186,191,212,230,232 ${wd}/DEST/utils/pnps.PrivateSNPs.sh
### sacct -j 17187383

module load intel/18.0 intelmpi/18.0 R/3.6.3

wd="/scratch/aob2x/dest"

#SLURM_ARRAY_TASK_ID=1


Rscript ${wd}/DEST/utils/pnps.privateSNPs.R ${SLURM_ARRAY_TASK_ID}
