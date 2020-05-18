#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:20:00 ### 0.5 hours
##SBATCH --mem 3G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.%A_%a.err # Standard error
#SBATCH -p biol8083
#SBATCH --account berglandlab

module load parallel

wd="/scratch/aob2x/dest"

moveFile () {

  #pop=B
  #chr=2L

  pop=${1}
  wd="/scratch/aob2x/dest"


  [ ! -d /project/berglandlab/DEST/dest_mapped/${pop}/ ] && mkdir /project/berglandlab/DEST/dest_mapped/${pop}

  rsync ${wd}/dest/wholeGenomeSyncData/${pop}_Chr*.gSYNC.gz* /project/berglandlab/DEST/dest_mapped/${pop}/

  echo $pop

}
export -f moveFile

parallel -j 1 moveFile ::: $( cat ${wd}/DEST/populationInfo/dpgp.ind.use.csv | sed '1d' | cut -f1 -d',' | sort | uniq )
