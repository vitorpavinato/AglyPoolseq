#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 0.5 hours
##SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

### sbatch ${wd}/DEST/misc_obsolete/move.10Nov2020.sh PoolSeq SNAPE NA NA 10Nov2020
# sacct -j 17962317


#module load htslib bcftools intel/18.0 intelmpi/18.0 parallel R/3.6.3


popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
#maf=NA; mac=NA; popSet="PoolSeq"; method="SNAPE"; version="10Nov2020"

wd="/scratch/aob2x/dest"

#echo "tabix"
#tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.ann.vcf.gz
#

echo "rsync"
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.gds /project/berglandlab/DEST/gds/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.gds
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz /project/berglandlab/DEST/vcf/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz.tbi /project/berglandlab/DEST/vcf/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz.tbi
rsync ${wd}/dest.${popSet}.${method}.${maf}.${mac}.header.bcf /project/berglandlab/DEST/bcf/dest.${popSet}.${method}.${maf}.${mac}.${version}.header.bcf
