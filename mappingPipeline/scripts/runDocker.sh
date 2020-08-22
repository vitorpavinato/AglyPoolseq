#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
##SBATCH --ntasks-per-node=10 # one core
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### most jobs should run in 60 minutes or less; the mitochondria takes a lot longer to run through pool-snp
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/dockerMap.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/dockerMap.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### test run as: sbatch --array=10 ${wd}/DEST/mappingPipeline/scripts/runDocker.sh
# sacct -j 14446076
# cat /scratch/aob2x/dest/slurmOutput/dockerMap.14446076_10.out
### modules
  module load singularity

### define a few things
  outputDir=/project/berglandlab/DEST/dest_mapped/pipeline_output
  wd=/scratch/aob2x/dest

### check
  #cat ${wd}/DEST/populationInfo/samps.csv | cut -f1,14 -d',' | grep -v "NA" | wc -l
  #ls -lh /project/berglandlab/DEST/dest_mapped/pipeline_output | wc -l

### get job number
  #SLURM_ARRAY_TASK_ID=10
  pop=$( cat ${wd}/DEST/populationInfo/samps.csv | cut -f1,14 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  srx=$( cat ${wd}/DEST/populationInfo/samps.csv | cut -f1,14 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )

  echo $pop
  echo $srx

  touch /scratch/aob2x/fastq/${srx}_1.fastq.gz
  touch /scratch/aob2x/fastq/${srx}_2.fastq.gz

### run docker
  singularity run \
  ${wd}/dmelsync_hpc.sif \
  /scratch/aob2x/fastq/${srx}_1.fastq.gz \
  /scratch/aob2x/fastq/${srx}_2.fastq.gz \
  ${pop} \
  ${outputDir} \
  --cores $SLURM_CPUS_PER_TASK \
  --max-cov 0.95 \
  --min-cov 4 \
  --base-quality-threshold 25 \
  --dont-prep \
  --do-snape \
  --do_poolsnp \
