#!/usr/bin/env bash
#
#SBATCH -J lme4qtl # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=10 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 7:00:00 # Running time of 1 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch /scratch/aob2x/dest/DEST/utils/get_geoEndemic.sh
# sacct -j


zcat /scratch/aob2x/dest/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | grep -v "#" | awk '{
npoly=0
nmissing=0
pops=""
if(length($5)==1) {
  for(i=10; i<=NF; i++) {
    split($i, sp, ":")
    if(sp[1]=="0/1") npoly++
    if(sp[1]=="0/1") pops=pops";"i-9
    if(sp[1]=="./.") nmissing++
    }
     print $1"\t"$2"\t"npoly"\t"nmissing"\t"$4"\t"$5"\t"pops
  }
}' > /scratch/aob2x/dest/geo_endemic/geo_endemic.delim
