#!/usr/bin/env bash
#SBATCH --job-name=meanCov # A single job name for the array
##SBATCH --ntasks-per-node=10 # one core
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00
#SBATCH --mem 4G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/meanCov.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/meanCov.%A_%a.err # Standard error
#SBATCH --account PAS1715

### modules
  module load samtools/1.10

### define a few things
  outputDir=/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output
  wd=/fs/scratch/PAS1715/aphidpool

### get job number
  pop=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  
  echo $pop
  
  touch /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/${pop}/${pop}.agly.bam
  
### run weGot_thisCoveraged.sh
  
  sh ${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/weGot_thisCoveraged.sh \
  /fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/${pop}/${pop}.agly.bam \
  ${pop} \
  ${outputDir}
