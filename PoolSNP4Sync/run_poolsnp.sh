#!/usr/bin/env bash
#
#SBATCH -J split_and_run # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### get bam files

### split on chromosome chunks

module load htslib bcftools parallel


## working & temp directory
  wd="/scratch/aob2x/dest"
  syncPath1="/project/berglandlab/DEST/dest_mapped/*/*masked.sync.gz"
  syncPath2="/project/berglandlab/DEST/dest_mapped/*/*/*masked.sync.gz"

## get job
  #SLURM_ARRAY_TASK_ID=4
  job=$( cat ${wd}/dest/poolSNP_jobs.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" )

  echo $job

## set up RAM disk
  tmpdir="/scratch/aob2x/test"

## get sub section
  subsection () {
    syncFile=${1}
    job=${2}
    jobid=$( echo ${job} | sed 's/,/_/g' )
    tmpdir=${3}

    pop=$( echo ${syncFile} | rev | cut -f1 -d'/' | rev | cut -f1 -d '.' )


    #syncFile=/project/berglandlab/DEST/dest_mapped/B/B.masked.sync.gz
    #job=2L,1,13766

    chr=$( echo $job | cut -f1 -d',' )
    start=$( echo $job | cut -f2 -d',' )
    stop=$( echo $job | cut -f3 -d',' )

    echo ${pop}_${jobid}

    tabix -b 2 -s 1 -e 2 \
    ${syncFile} \
    ${chr}:${start}-${stop} > ${tmpdir}/${pop}_${jobid}

  }
  export -f subsection

  parallel -j 1 subsection ::: $( ls ${syncPath1} ${syncPath2} ) ::: ${job} ::: ${tmpdir}

### paste function (done with R)
  Rscript --no-save --no-restore paste.R ${job} ${tmpdir}
