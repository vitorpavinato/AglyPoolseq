#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 0.5 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/index.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/index.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch --array=1-6 /scratch/aob2x/dest/DEST/misc_obsolete/indexBamFiles.sh
# sacct -j 13715912
module load samtools

melBamPath="/project/berglandlab/DEST/dest_mapped/0803_three_redos/*/*mel.bam"
simBamPath="/project/berglandlab/DEST/dest_mapped/0803_three_redos/*/*sim.bam"

bams=$( ls ${melBamPath} ${simBamPath} )

#SLURM_ARRAY_TASK_ID=1

# echo ${bams} | tr ' ' '\n' | wc -l ; shoudl be 492

bamUse=$( echo ${bams} | tr ' ' '\n' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

samtools index $bamUse
