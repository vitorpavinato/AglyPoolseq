#!/usr/bin/env bash
#SBATCH --job-name=meanCovDepthR # A single job name for the array
##SBATCH --ntasks-per-node=10 # one core
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00
#SBATCH --mem 4G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/meanCovDepthR.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/meanCovDepthR.%A_%a.err # Standard error
#SBATCH --account PAS1715

### modules
  module load R/4.0.2-gnu9.1 


### define a few things
  outputDir=/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/aggregated
  wd=/fs/scratch/PAS1715/aphidpool

### get job number
  pop=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_aggregated.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  
  echo $pop
  
### run mean_coverage_across_runs.R  
  Rscript --vanilla ${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/mean_coverage_across_runs.R \
  ${outputDir} \
  ${pop} \
  ${outputDir}
