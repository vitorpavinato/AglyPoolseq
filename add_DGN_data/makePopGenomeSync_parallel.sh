#!/usr/bin/env bash
#
#SBATCH -J makePopSync # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 20 # on one node
#SBATCH -t 4:30:00 ### 5 minutes
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/makePopSync.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/makePopSync.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load htslib

wd="/scratch/aob2x/dest"

#####################
### get jobs task ###
#####################

#SLURM_ARRAY_TASK_ID=1
##chr_i=$( echo "${SLURM_ARRAY_TASK_ID}%5+1" | bc )
##pop_i=$( echo "${SLURM_ARRAY_TASK_ID}%35+1" | bc )

pop=$( grep  "^${SLURM_ARRAY_TASK_ID}[[:space:]]" /scratch/aob2x/dest/dgn/pops.delim | cut -f3 )
chr=$( grep  "^${SLURM_ARRAY_TASK_ID}[[:space:]]" /scratch/aob2x/dest/dgn/pops.delim | cut -f2 )


echo $pop
echo $chr

############################
### paste per chromosome ###
############################

#if [ ! -f /scratch/aob2x/dest/dgn/csvData/${pop}_Chr${chr}.csv ]; then
#
#    paste -d',' /scratch/aob2x/dest/dgn/longData/${pop}_*_Chr${chr}*long  \
#    /scratch/aob2x/dest/dgn/csvData/${pop}_Chr${chr}.csv
#
#fi

## pop=RAL; chr=2L

samps=$( grep "^"${pop}"," < ${wd}/DEST/populationInfo/dpgp.ind.use.csv | cut -f2 -d',' | sed "s/$/_Chr${chr}/g" | tr '\n' '|' | sed 's/.$//' )

nFiles_obs=$( ls /scratch/aob2x/dest/dgn/longData/${pop}_*_Chr${chr}*long | grep -E "${samps}" | wc -l )

nFiles_exp=$( grep "^"${pop}"," < ${wd}/DEST/populationInfo/dpgp.ind.use.csv | cut -f9 -d',' | head -n1 )

echo ${pop}","${chr}","${nFiles_obs}","${nFiles_exp} >> /scratch/aob2x/dest/dgn/confirm_files



### use tmpdir
#SLURM_ARRAY_TASK_ID=1; SLURM_JOB_ID=1; SLURM_NTASKS_PER_NODE=1
tmpdir=/dev/shm/$USER/
[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}


### first make sub_input files
ls /scratch/aob2x/dest/dgn/longData/${pop}_*_Chr${chr}*long | grep -E "${samps}" | split -l 20 - /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/jobs


subPaste () {
  job=${1} # job=jobsaa
  tmpdir=${2} # tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  paste -d',' $( cat ${tmpdir}/${job} ) | \
  awk -F' ' -v chr=${chr} '
  {
  nN=gsub(/N/,"",$2)
  nA=gsub(/A/,"",$2)
  nT=gsub(/T/,"",$2)
  nC=gsub(/C/,"",$2)
  nG=gsub(/G/,"",$2)

  nObs=nA+nT+nC+nG

  print nA":"nT":"nC":"nG":"nN":0"

}' > ${tmpdir}/sync.${job}.sync

}
export -f subPaste

cd /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}
parallel -j ${SLURM_NTASKS_PER_NODE} subPaste ::: $( ls jobs* | grep -v "csv" ) ::: ${tmpdir}


### add sync
paste $( ls sync*sync ) | awk '{
  nA=0
  nT=0
  nC=0
  nG=0
  nN=0


  for(i=1; i<=NF; i++) {
    spN=split($i, sp, ":")
    nA=nA+sp[1]
    nT=nT+sp[2]
    nC=nC+sp[3]
    nG=nG+sp[4]
    nN=nN+sp[5]
  }

  print nA":"nT":"nC":"nG":"nN":0"


}' > tmp.gSYNC


paste -d' ' \
/scratch/aob2x/dest/referenceGenome/r5/${chr}.long \
tmp.gSYNC | \
awk -F' ' -v chr=${chr} '
{

print chr"\t"NR"\t"toupper($1)"\t"$2

}' | bgzip -c > /scratch/aob2x/dest/dest/wholeGenomeSyncData/${pop}_Chr${chr}.sync.gz


rm -fr /dev/shm/$USER/${SLURM_JOB_ID}/*

#tabix -f -b 2 -s 1 -e 2 /scratch/aob2x/dest/dest/wholeGenomeSyncData/${pop}_Chr${chr}.sync.gz
