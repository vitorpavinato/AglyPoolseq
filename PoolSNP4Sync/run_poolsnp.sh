#!/usr/bin/env bash
#
#SBATCH -J download_DGN # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/download_DGN.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/download_DGN.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### get bam files

### split on chromosome chunks

module load htslib bcftools


## working & temp directory
  wd="/scratch/aob2x/dest"
  tmpdir="/scratch/aob2x/test"

## get job
  #SLURM_ARRAY_TASK_ID=4
  job=$( cat ${wd}/dest/poolSNP_jobs.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" )

  echo $job

## get sub section
  subsection () {
    syncFile=${1}
    job=${2}

    #syncFile=/project/berglandlab/DEST/dest_mapped/B/B.masked.sync.gz
    #job=2L,1,13766

    chr=$( echo $job | cut -f1 -d',' )
    start=$( echo $job | cut -f2 -d',' )
    stop=$( echo $job | cut -f3 -d',' )

    tabix -b 2 -s 1 -e 2 \
    ${syncFile} \
    ${chr}:${start}-${stop}

  }
