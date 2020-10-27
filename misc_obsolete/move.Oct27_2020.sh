#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 0.5 hours
##SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.16425773.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

### sbatch ${wd}/DEST/misc_obsolete/move.Oct27_2020.sh SNAPE 001 50
# sacct -j 17036246


module load htslib bcftools intel/18.0 intelmpi/18.0 parallel R/3.6.3


popSet=${1}
method=${2}
maf=${3}
mac=${4}
#maf=001; mac=50; popSet="PoolSeq"; method="SNAPE"

wd=/scratch/aob2x/dest

tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.ann.vcf.gz

wd="/scratch/aob2x/dest"
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.ann.gds /project/berglandlab/DEST/gds/.
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.ann.vcf.gz /project/berglandlab/DEST/vcf/.
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.ann.vcf.gz.tbi /project/berglandlab/DEST/vcf/.
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.header.bcf /project/berglandlab/DEST/bcf/.
