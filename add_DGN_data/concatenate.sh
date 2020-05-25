#!/usr/bin/env bash
#
#SBATCH -J download_DGN # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/concatenate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/concatenate.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1 



module load htslib bcftools parallel

### note: this script uses the ram disk to write temporary files.

wd="/scratch/aob2x/dest"

#SLURM_ARRAY_TASK_ID=10
pop=$( cat ${wd}/dgn/pops.delim | cut -f3 | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d" )



### un-compress

uncomp () {
  echo uncompressing: ${1}

  outfn=$( echo ${1} | sed 's/.gz//g' )
  bgzip -c -d ${1} > ${outfn}
}
export -f uncomp

parallel uncomp ::: $( ls ${wd}/dest/wholeGenomeSyncData/${pop}_*.gSYNC.gz )

### change out missing data annotation
swapfun() {
  echo swapping: ${1}
  sed -i 's/0:0:0:0:0:0/.:.:.:.:.:./g' ${1}
}
export -f swapfun

parallel swapfun ::: $( ls ${wd}/dest/wholeGenomeSyncData/${pop}_*.gSYNC )

### concatenate
cat \
${wd}/dest/wholeGenomeSyncData/${pop}_2L.gSYNC \
${wd}/dest/wholeGenomeSyncData/${pop}_2R.gSYNC \
${wd}/dest/wholeGenomeSyncData/${pop}_3L.gSYNC \
${wd}/dest/wholeGenomeSyncData/${pop}_3R.gSYNC \
${wd}/dest/wholeGenomeSyncData/${pop}_X.gSYNC |
bgzip -c >
${wd}/dest/wholeGenomeSyncData/${pop}_X.sync.gz
