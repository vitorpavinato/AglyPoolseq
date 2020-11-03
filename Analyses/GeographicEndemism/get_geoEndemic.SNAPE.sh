#!/usr/bin/env bash
#
#SBATCH -J geoEndemic # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 7:00:00 # Running time of 1 hours
#SBATCH --mem 1G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch --array=1-12 /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/get_geoEndemic.SNAPE.sh
# sacct -j 18548987

# SLURM_ARRAY_TASK_ID=5

maf=$( echo "0, 0.001, 0.003, 0.005, 0.007, 0.01, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3" | sed 's/ //g' | cut -d',' -f ${SLURM_ARRAY_TASK_ID} )

zcat /scratch/aob2x/dest/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | grep -v "#" | awk -v maf=${maf} '{
npoly=0
nmissing=0
af=0
pops=""
if(length($5)==1) {
  for(i=10; i<=NF; i++) {
    split($i, sp, ":")
    if(sp[1]=="0/1" && sp[5]>maf && sp[5]<(1-maf)) npoly++
    if(sp[1]=="0/1" && sp[5]>maf && sp[5]<(1-maf)) pops=pops";"i-9
    if(sp[1]=="0/1" && sp[5]>maf && sp[5]<(1-maf)) af=af"+"sp[5]
    if(sp[1]=="./.") nmissing++
    }
    if(npoly>1) print maf"\t"$1"\t"$2"\t"npoly"\t"nmissing"\t"$4"\t"$5"\t"pops"\t"af
  }
}' | head > /scratch/aob2x/dest/geo_endemic/geo_endemic.SNAPE.${SLURM_ARRAY_TASK_ID}.delim
