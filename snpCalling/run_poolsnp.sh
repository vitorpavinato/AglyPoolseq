#!/usr/bin/env bash
#
#SBATCH -J split_and_run # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00 ### 15 minutes
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch --array=1-$( wc -l ${wd}/dest/poolSNP_jobs.csv | cut -f1 -d' ' ) ${wd}/DEST/PoolSNP4Sync/run_poolsnp.sh
### sbatch --array=3000-3500 ${wd}/DEST/snpCalling/run_poolsnp.sh
### sacct -j 12761688

module load htslib bcftools parallel intel/18.0 intelmpi/18.0 R/3.6.0


## working & temp directory
  wd="/scratch/aob2x/dest"
  syncPath1="/project/berglandlab/DEST/dest_mapped/*/*masked.sync.gz"
  syncPath2="/project/berglandlab/DEST/dest_mapped/*/*/*masked.sync.gz"
  outdir="/scratch/aob2x/dest/sub_vcfs"

## get job
  #SLURM_ARRAY_TASK_ID=4
  job=$( cat ${wd}/dest/poolSNP_jobs.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" )
  jobid=$( echo ${job} | sed 's/,/_/g' )
  echo $job

## set up RAM disk
  ## rm /scratch/aob2x/test/*
  #tmpdir="/scratch/aob2x/test"
  #SLURM_JOB_ID=1; SLURM_ARRAY_TASK_ID=4
  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

## get sub section
  subsection () {
    syncFile=${1}
    job=${2}
    jobid=$( echo ${job} | sed 's/,/_/g' )
    tmpdir=${3}

    pop=$( echo ${syncFile} | rev | cut -f1 -d'/' | rev | sed 's/.masked.sync.gz//g' )


    #syncFile=/project/berglandlab/DEST/dest_mapped/pipeline_output/ES_ba_12_spring/ES_ba_12_spring.masked.sync.gz

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

  echo "subset"

  parallel -j 1 subsection ::: $( ls ${syncPath1} ${syncPath2} | grep -v "SNAPE" ) ::: ${job} ::: ${tmpdir}

### paste function
  echo "paste"
  Rscript --no-save --no-restore ${wd}/DEST/snpCalling/paste.R ${job} ${tmpdir}

### run through PoolSNP
  echo "poolsnp"
  cat ${tmpdir}/allpops.sites | python ${wd}/DEST/snpCalling/PoolSnp.py \
  --sync - \
  --min-cov 4 \
  --max-cov 0.95 \
  --min-count 4 \
  --min-freq 0.001 \
  --miss-frac 0.5 \
  --names $( cat ${tmpdir}/allpops.names | tr '\n' ',' | sed 's/,$//g' ) > ${tmpdir}/${jobid}.vcf

### compress and clean up
  echo "compress and clean"
  bgzip -c ${tmpdir}/${jobid}.vcf > ${outdir}/${jobid}.vcf.gz
  tabix -p vcf ${outdir}/${jobid}.vcf.gz

  #echo "vcf -> bcf "
  #bcftools view -Ou ${tmpdir}/${jobid}.vcf.gz > ${outdir}/${jobid}.bcf

  rm -fr ${tmpdir}

### done
  echo "done"
