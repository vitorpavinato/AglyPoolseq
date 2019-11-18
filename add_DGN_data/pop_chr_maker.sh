#!/usr/bin/env bash
#
#SBATCH -J pop_chr_maker # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:30:00 ### 30 minutes
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/pop_chr_maker.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/pop_chr_maker.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


ls /scratch/aob2x/dest/dgn/longData/* | rev | cut -f1 -d'/' | rev | cut -f1 -d'_' | sort | uniq | \
awk '{print "2L\t"$0"\n2R\t"$0"\n3L\t"$0"\n3R\t"$0"\nX\t"$0}' | awk '{print NR"\t"$0}' > \
/scratch/aob2x/dest/dgn/pops.delim
