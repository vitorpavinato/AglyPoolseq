#!/usr/bin/env bash
#
#SBATCH --job-name=download_SRA # A single job name for the array
##SBATCH --nodes=1 --ntasks-per-node=10 # one core
#SBATCH -c 2
#SBATCH -N 1 # on one node
#SBATCH --time=6:00:00 ### 6 hours
#SBATCH --mem 8G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/download_SRA.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/download_SRA.%A_%a.err # Standard error
#SBATCH --account PAS1715


### grep -E "ES_ba_12|AT_gr_12" /scratch/aob2x/dest/DEST/populationInfo/samps.csv | cut -f1,13 -d',' > /scratch/aob2x/fastq/todl.csv
### run as: sbatch --array=1-$( wc -l < /scratch/aob2x/fastq/todl.csv ) /scratch/aob2x/dest/DEST/mappingPipeline/scripts/download_SRA.sh
### sacct -j 18750292
module load sratoolkit/2.10.7

wd=/fs/scratch/PAS1715/aphidpool

#SLURM_ARRAY_TASK_ID=1

sranum=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_combined.csv | cut -f1,14 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )
sampName=$( cat ${wd}/DEST-AglyPoolseq/populationInfo/fieldPools_combined.csv | cut -f1,14 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )

echo $sampName " / " $sranum

if [ ! -d ${wd}/fastq ]; then
  mkdir ${wd}/fastq
fi

prefetch \
-o /fs/scratch/PAS1715/aphidpool/fastq/${sranum}.sra \
${sranum}

fasterq-dump \
--split-files \
--split-3 \
--outfile /fs/scratch/PAS1715/aphidpool/fastq/${sranum} \
-e 10 \
/fs/scratch/PAS1715/aphidpool/fastq/${sranum}.sra

gzip /fs/scratch/PAS1715/aphidpool/fastq/${sranum}_1.fastq
gzip /fs/scratch/PAS1715/aphidpool/fastq/${sranum}_2.fastq

#rm /fs/scratch/PAS1715/aphidpool/fastq/${sranum}.sra
