#!/usr/bin/env bash
#
#SBATCH -J geoEndemic # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 00:30:00 # Running time of 1 hours
#SBATCH --mem 60G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab

# sbatch /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.collect.sh
# sacct -j 18781752
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.18781752
module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3


# SLURM_ARRAY_TASK_ID=9
#maf=$( cat /scratch/aob2x/dest/geo_endemic/jobs.txt | cut -f1 | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d" )


#Rscript /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.collect.goodSamps.SNAPE.R
Rscript /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/geographic_clustering_endemic.collect.PoolSeq.PoolSNP.R
