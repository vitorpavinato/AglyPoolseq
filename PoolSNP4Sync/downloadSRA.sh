#!/bin/bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 ### 36 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/download_SRA.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/download_SRA.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# run as: sbatch --array=1-$( wc -l < ${wd}/DEST/populationInfo/samps.csv ) ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
# demo as: sbatch --array=1-5 ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh


module load sratoolkit

### set working directory
cd /scratch/aob2x/dest/

### define a few paths
wd=/scratch/aob2x/dest
fastq_directory=/scratch/aob2x/dest/fastq

# SLURM_ARRAY_TASK_ID=175

### get SRA accession number
sra_col=$( cat ${wd}/DEST/populationInfo/samps.csv | head -n1 | tr ',' '\n' | awk '{print NR"\t"$0}' | grep "SRA_accession" | cut -f1 )
sra_num=$( cat ${wd}/DEST/populationInfo/samps.csv | awk -v samp=${SLURM_ARRAY_TASK_ID} -v sraCol=${sra_col} -F, '{if((NR-1)==samp) print $sraCol}' )


if [ ${sra_num}!=NA ]; then
  fastq-dump \
  --split-files \
  -O ${fastq_directory} \
  --gzip \
  -Q 33 \
  --disable-multithreading \
  ${sra_num}
fi
