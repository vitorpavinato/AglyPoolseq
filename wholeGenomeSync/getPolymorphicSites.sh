#!/usr/bin/env bash
#
#SBATCH -J getPolymorphicSites # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00 ### 15 minutes
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/getPolymorphicSites.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/getPolymorphicSites.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load htslib
module load gcc/7.1.0 bedtools/2.26.0

### pop=UM; chr=2L

pop=${1}
chr=${2}
in_dir=/scratch/aob2x/dest/dest/wholeGenomeSyncData
out_dir=/scratch/aob2x/dest/dest/sites

bgzip -cd ${in_dir}/${pop}_Chr${chr}.sync.gz  |
awk -v MAC=1 -v pop=${pop} -v chr=${chr} -v od=${out_dir} '{

  split($4, alleleCounts, ":")

  numAlleles=0
  for(i=1; i<=4; i++) {
    if(alleleCounts[i]>=MAC) {
      numAlleles++
    }
  }
  print $1"\t"$2"\t"$2+1"\t"numAlleles > od/pop

}'
