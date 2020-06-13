#!/usr/bin/env bash
#
#SBATCH -J touchup # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 48:00:00 ### 6 hours
#SBATCH --mem 30G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/touchup.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/touchup.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

wd=/scratch/aob2x/dest
### run as: sbatch --array=1-$( tail -n1 /scratch/aob2x/fastq/todl.csv | cut -f3 -d',' ) ${wd}/DEST/misc_obsolete/runTouchup.sh

module load singularity

#SLURM_ARRAY_TASK_ID=1

sranum=$( grep ",${SLURM_ARRAY_TASK_ID}$" /scratch/aob2x/fastq/todl.csv | cut -f2 -d',' )
sampName=$( grep ",${SLURM_ARRAY_TASK_ID}$" /scratch/aob2x/fastq/todl.csv | cut -f1 -d',' )

echo $sampName " / " $sranum


cd /scratch/aob2x/test

singularity run \
/scratch/aob2x/dest/dmelsync_hpc.sif \
/scratch/aob2x/fastq/${sranum}_1.fastq.gz \
/scratch/aob2x/fastq/${sranum}_1.fastq.gz \
${sampName}_sub \
${sampName} \
-c 20
