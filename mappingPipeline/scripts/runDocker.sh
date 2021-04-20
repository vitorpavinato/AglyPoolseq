#!/usr/bin/env bash
#SBATCH --job-name=dockerMap # A single job name for the array
##SBATCH --ntasks-per-node=10 # one core
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 03:30:00 ### most jobs should run in 60 minutes or less; the mitochondria takes a lot longer to run through pool-snp
#SBATCH --mem 12G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/dockerMap.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/dockerMap.%A_%a.err # Standard error
#SBATCH --account PAS1715

### test run as: sbatch --array=10 ${wd}/DEST-AglyPoolseq/mappingPipeline/scripts/runDocker.sh
# sacct -j 14446076
# cat /fs/scratch/PAS1715/aphidpool/slurmOutput/dockerMap.16904399_2.out | grep -v "Reference is N; most frequent allele is calculated in position" | less -S
# cat /fs/scratch/PAS1715/aphidpool/slurmOutput/dockerMap.16904399_2.out

### modules
  module load singularity

### define a few things
  outputDir=/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output
  wd=/fs/scratch/PAS1715/aphidpool

### check
  #cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | wc -l
  #ls -lh ${wd}/dest_mapped/pipeline_output | wc -l

### get job number
  #SLURM_ARRAY_TASK_ID=10
  # for the technical replicated run files
  pop=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  prx=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )
  numFlies=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,12 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )

  echo ${SLURM_ARRAY_TASK_ID}
  echo $pop
  echo $prx
  echo $numFlies
  
  # for the replicated run files
  touch /fs/scratch/PAS1715/aphidpool/fastq/${prx}.R1.fastq.gz
  touch /fs/scratch/PAS1715/aphidpool/fastq/${prx}.R2.fastq.gz
  
  singularity run \
  ${wd}/aglypoolseq_latest.sif \
  /fs/scratch/PAS1715/aphidpool/fastq/${prx}_1.fastq.gz \
  /fs/scratch/PAS1715/aphidpool/fastq/${prx}_2.fastq.gz \
  ${pop} \
  ${outputDir} \
  --cores $SLURM_CPUS_PER_TASK \
  --min-indel 5 \
  --max-cov 0.95 \
  --min-cov 3 \
  --base-quality-threshold 25 \
  --num-flies ${numFlies} \
  --dont-prep \
  --do_poolsnp
