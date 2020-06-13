#!/usr/bin/env bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load sratoolkit/2.9.1

wd=/scratch/aob2x/dest
#SLURM_ARRAY_TASK_ID=225

sranum=$( sed -n "${SLURM_ARRAY_TASK_ID}p" ${wd}/DEST/populationInfo/samps.csv | cut -f13 -d',' )

prefetch \
-O /scratch/aob2x/${sranum} \
${sranum}

fasterq-dump \
--split-files \
--outfile /scratch/aob2x/temp_${sranum} \
-e 10 \
/scratch/aob2x/${sranum}.sra

gzip /scratch/aob2x/dest/fastq/${sranum}_1.fastq
gzip /scratch/aob2x/dest/fastq/${sranum}_2.fastq
