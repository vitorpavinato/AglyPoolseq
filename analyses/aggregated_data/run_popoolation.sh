#!/usr/bin/env bash
#SBATCH -J run_popoolations
#SBATCH -c 2
#SBATCH -N 1 # on one node
#SBATCH -t 02:00:00 
#SBATCH --mem 8G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/run_popoolation.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/run_popoolation.%A_%a.err # Standard error
#SBATCH --account PAS1715

## Load modules
module load samtools/1.9

## Set the running environment
# PATHS
wd="/fs/scratch/PAS1715/aphidpool"
inputDir=/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/aggregated
outputDir=/fs/scratch/PAS1715/aphidpool/popoolation_sumstats

# POOL NAMES AND SIZE
poolname=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
poolcode=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )
poolsize=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,12 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )

## Run monitoring outputs
echo ${SLURM_ARRAY_TASK_ID}
echo $poolname
echo $poolcode
echo $poolsize

## Run Popoolation summary statistics for each pool
echo "running poopolation summary statistics"
sh ${wd}/DEST-AglyPoolseq/analyses/aggregated_data/popoolation_summStats.sh ${inputDir} ${outputDir} ${poolname} ${poolsize}

echo "done"
