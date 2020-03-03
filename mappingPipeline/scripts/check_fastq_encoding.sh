#!/bin/bash
#
#SBATCH -J check_fastq_encoding # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/check_fastq_encoding.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/check_fastq_encoding.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


gunzip SRX2885350_2.fastq.gz | head awk 'NR % 4 == 0' | guess-encoding.py -n 1000
