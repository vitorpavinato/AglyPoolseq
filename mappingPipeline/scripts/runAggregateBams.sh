#!/usr/bin/env bash
#SBATCH --job-name=aggregateBams # A single job name for the array
##SBATCH --ntasks-per-node=10 # one core
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 02:00:00
#SBATCH --mem 20G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/aggregateBams.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/aggregateBams.%A_%a.err # Standard error
#SBATCH --account PAS1715

### modules
  #module load singularity
  module load samtools/1.9

### define a few things
  wd=/fs/scratch/PAS1715/aphidpool
  inputDir=/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output
  outputDir=/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/aggregated

### get job number
  pool=$( cat ${wd}/AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  
  echo $pool
  
  touch /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/${pool}_*/${pool}_*.sorted_merged.bam
  
### run aggregateBams.sh
  
  sh ${wd}/AglyPoolseq/mappingPipeline/scripts/aggregate_bams.sh \
  sorted_merged \
  ${inputDir} \
  ${pool} \
  ${outputDir} \
  --cores $SLURM_CPUS_PER_TASK
