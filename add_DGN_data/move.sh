#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:20:00 ### 0.5 hours
##SBATCH --mem 3G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

module load parallel htslib

wd="/scratch/aob2x/dest"

# rm -R -v !("/project/berglandlab/DEST/dest_mapped/example_pipeline_output")


index_moveFile () {

  #pop=B
  #chr=2L

  pop=${1}
  wd="/scratch/aob2x/dest"


  [ ! -d /project/berglandlab/DEST/dest_mapped/${pop}/ ] && mkdir /project/berglandlab/DEST/dest_mapped/${pop}

  tabix -f -b 2 -s 1 -e 2 ${wd}/dest/wholeGenomeSyncData/${pop}.sync.gz

  rsync ${wd}/dest/wholeGenomeSyncData/${pop}.sync.gz* /project/berglandlab/DEST/dest_mapped/${pop}/

  echo $pop

}
export -f index_moveFile

parallel -j 1 index_moveFile ::: $( cat ${wd}/DEST/populationInfo/dpgp.ind.use.csv | sed '1d' | cut -f1 -d',' | sort | uniq )
