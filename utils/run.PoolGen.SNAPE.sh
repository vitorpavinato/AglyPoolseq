#!/usr/bin/env bash
#
#SBATCH -J poolgen # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 1 hours
#SBATCH --mem 150G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/poolgen.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/poolgen.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account berglandlab


### run as: sbatch --array=1-246 ${wd}/DEST/utils/run.PoolGen.SNAPE.sh
### sacct -j 19362455
### cat /scratch/aob2x/dest/slurmOutput/poolgen.19358567_1.out

### cat /scratch/aob2x/dest/slurmOutput/poolgen.19358314_1.out


module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 parallel bcftools

export wd="/scratch/aob2x/dest"

### get sample information
  #SLURM_ARRAY_TASK_ID=1
  export popName=$( cat ${wd}/DEST/populationInfo/samps.csv | grep "pooled" | cut -f1,14 -d',' | grep -v "NA"  | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  export numFlies=$( cat ${wd}/DEST/populationInfo/samps.csv | grep "pooled" | cut -f1,12 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )


## set up RAM disk
echo "make ram disk"
  ## rm /scratch/aob2x/test/*
  #tmpdir="/scratch/aob2x/test"
  #SLURM_JOB_ID=1
  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  export tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

### run each chromosome separately
echo "run"

  runPoolGen () {
    chr=${1}
    ###chr=2L
    #tabix -p vcf /project/berglandlab/DEST/vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz ${chr} | head -n 100000 > ${tmpdir}/${chr}.vcf
    echo ${chr} ${popName} ${numFlies}

    echo "extracting bed"
      zcat -c /project/berglandlab/DEST/dest_mapped/pipeline_output/${popName}/${popName}.SNAPE.complete.bed.gz | grep "${chr}" > \
      ${tmpdir}/${popName}.${chr}.bed

    echo "extracting chromosome"

      bcftools view \
      -O v \
      /project/berglandlab/DEST/vcf/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz ${chr} > ${tmpdir}/${chr}.vcf

    if [[ $chr == "X" ]]; then
      numChr=${numFlies}
      echo ${numChr}
    else
      #numChr =
      numChr=$(( $numFlies*2 ))
      echo ${numChr}
    fi
    echo "numChr:" ${numChr}

    echo "running PoolGen"
      python3 ${wd}/DEST/utils/PoolGen.py \
      --input ${tmpdir}/${chr}.vcf \
      --step 200000 \
      --window 200000 \
      --pool-size ${numChr} \
      --min-sites-frac 0.75 \
      --BED ${tmpdir}/${popName}.${chr}.bed \
      --sample ${popName} \
      --output /scratch/aob2x/dest/PoolGenOut_SNAPE/${popName}.${chr}

    rm ${tmpdir}/${chr}.vcf

  }
  export -f runPoolGen

  parallel -j 5 runPoolGen ::: 2L 2R 3L 3R X


### clean up
  rm -fr ${tmpdir}
