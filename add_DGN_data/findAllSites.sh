#!/usr/bin/env bash
#
#SBATCH -J findAllSites # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/findAllSites.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/findAllSites.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


union_bed () {
  chr=${1}
  dir=${2}

  cat ${dir}/*${chr}*.sync | cut -f2 | sort | uniq | sort -k1,1n | awk -v chr=${chr} '{print "chr"chr"\t"$0"\t"$0+1}' > /mnt/spicy_2/dest/african/dpgp_sites/dpgp_${chr}.bed
}
export -f union_bed


parallel --gnu -j1 union_bed ::: 2L 2R 3L 3R X ::: /mnt/spicy_2/dest/african/+(dpgp2|dpgp3)_sequences_temp
