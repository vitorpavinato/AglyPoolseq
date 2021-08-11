#!/usr/bin/env bash
#SBATCH -J split_and_run # A single job name for the array
#SBATCH -c 3
#SBATCH -N 1 # on one node
#SBATCH -t 3:00:00 ### most jobs should run in 60 minutes or less; the mitochondria takes a lot longer to run through pool-snp
#SBATCH --mem 12G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH --account PAS1715

## Load modules
module load htslib
module load bcftools/1.9.2
module load vcftools/0.1.16
module load parallel2/19.10
module load R/4.0.2-gnu9.1
module load python/3.6

## working & temp directory
  wd="/fs/scratch/PAS1715/aphidpool"
  #outdir="/fs/scratch/PAS1715/aphidpool/sub_vcfs"
  outdir="/fs/scratch/PAS1715/aphidpool/sub_vcfs_aggregated"
  popSet=${1}
  method=${2}
  maf=${3}
  mac=${4}
  version=${5}
  jobs=${6}
  #maf=01; mac=50; popSet="all"; method="PoolSNP"; version="paramTest"; jobs=poolSNP_jobs.sample.csv
  
## Check if outdir exists
  [ ! -d ${outdir} ] && mkdir ${outdir}

## get list of SNYC files based on popSet & method - all replicates of each pool vs pools of aggregated replicates
### full list
  syncPath1orig="/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/*/*masked.sync.gz"
  syncPath2orig="/fs/scratch/PAS1715/aphidpool/dest_mapped/pipeline_output/aggregated/*/*masked.sync.gz"
  
### target pops
  if [[ "${popSet}" == "PoolSeq" ]]; then
    syncPath1=""
    syncPath2=${syncPath2orig}
  elif [[ "${popSet}" == "all" ]]; then
    syncPath1=${syncPath1orig}
    syncPath2=""
  fi
  
## get job
  job=$( cat ${wd}/${jobs} | sed "${SLURM_ARRAY_TASK_ID}q;d" )
  jobid=$( echo ${job} | sed 's/,/_/g' )
  echo $job

# set up RAM disk
  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

# set up RAM disk
#  [ ! -d /fs/scratch/PAS1715/aphidpool/$USER/ ] && mkdir /fs/scratch/PAS1715/aphidpool/$USER/
#  [ ! -d /fs/scratch/PAS1715/aphidpool/$USER/${SLURM_JOB_ID} ] && mkdir /fs/scratch/PAS1715/aphidpool/$USER/${SLURM_JOB_ID}
#  [ ! -d /fs/scratch/PAS1715/aphidpool/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /fs/scratch/PAS1715/aphidpool/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}
#
#  tmpdir=/fs/scratch/PAS1715/aphidpool/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

## get sub section
  subsection () {
    syncFile=${1}
    job=${2}
    jobid=$( echo ${job} | sed 's/,/_/g' )
    tmpdir=${3}

    pop=$( echo ${syncFile} | rev | cut -f1 -d'/' | rev | sed 's/.masked.sync.gz//g' )
    
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
  
  if [[ "${method}" == "SNAPE" ]]; then
    echo "SNAPE" ${method}
    parallel -j 1 subsection ::: $( ls ${syncPath1} ${syncPath2} | grep "SNAPE" | grep "monomorphic" ) ::: ${job} ::: ${tmpdir}
  elif [[ "${method}" == "PoolSNP" ]]; then
    echo "PoolSNP" ${method}
    parallel -j 1 subsection ::: $( ls ${syncPath1} ${syncPath2} | grep -v "SNAPE" ) ::: ${job} ::: ${tmpdir}
  fi

### paste function
  echo "paste"
  Rscript --no-save --no-restore ${wd}/DEST-AglyPoolseq/snpCalling/paste.R ${job} ${tmpdir} ${method}

### run through PoolSNP
  echo "poolsnp"

  if [[ "${method}" == "SNAPE" ]]; then
    echo $method
    cat ${tmpdir}/allpops.${method}.sites | python ${wd}/DEST-AglyPoolseq/snpCalling/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.99 \
    --miss-frac 0.5 \
    --min-count 0 \
    --min-freq 0 \
    --posterior-prob 0.9 \
    --SNAPE \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf

  elif [[ "${method}"=="PoolSNP" ]]; then
    echo $method

    cat ${tmpdir}/allpops.${method}.sites | python ${wd}/DEST-AglyPoolseq/snpCalling/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.99 \
    --min-count ${mac} \
    --min-freq 0.${maf} \
    --miss-frac 0.5 \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  fi

### compress and clean up
  echo "compress and clean"

  cat ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf | vcf-sort | bgzip -c > ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  tabix -p vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

### done
  echo "done"
