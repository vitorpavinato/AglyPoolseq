#!/usr/bin/env bash
#
#SBATCH -J geoEndemic # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 2:00:00 # Running time of 1 hours
#SBATCH --mem 1G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch --array=1-50 /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/get_geoEndemic.SNAPE.sh
# sacct -j 18549212
# sbatch --array=38 /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/get_geoEndemic.SNAPE.sh
# sacct -j 18550136

# sbatch /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/get_geoEndemic.SNAPE.sh



# SLURM_ARRAY_TASK_ID=5
module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3 htslib bcftools

wd=/scratch/aob2x/dest

### fix MAF at 0.05
bcftools view \
-S /scratch/aob2x/dest/DEST/Analyses/GeographicEndemism/goodSamps.delim \
${wd}/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz | \
grep -v "#" | awk -v maf=0.05 '{
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
    if(npoly>0) print maf"\t"$1"\t"$2"\t"npoly"\t"nmissing"\t"$4"\t"$5"\t"pops"\t"af
  }
}' > /scratch/aob2x/dest/geo_endemic/geo_endemic.SNAPE.goodSamps.delim






### runs across multiple MAF thresholds
#Rscript ${wd}/DEST/Analyses/GeographicEndemism/makeJobs.R
#
#maf=$( cat ${wd}/geo_endemic/jobs.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1)
#chr=$( cat ${wd}/geo_endemic/jobs.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" | sed 's/"//g' | cut -f2)
#
#echo ${maf} ${chr}
#
#tabix ${wd}/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz ${chr} | grep -v "#" | awk -v maf=${maf} '{
#npoly=0
#nmissing=0
#af=0
#pops=""
#if(length($5)==1) {
#  for(i=10; i<=NF; i++) {
#    split($i, sp, ":")
#    if(sp[1]=="0/1" && sp[5]>maf && sp[5]<(1-maf)) npoly++
#    if(sp[1]=="0/1" && sp[5]>maf && sp[5]<(1-maf)) pops=pops";"i-9
#    if(sp[1]=="0/1" && sp[5]>maf && sp[5]<(1-maf)) af=af"+"sp[5]
#    if(sp[1]=="./.") nmissing++
#    }
#    if(npoly>1) print maf"\t"$1"\t"$2"\t"npoly"\t"nmissing"\t"$4"\t"$5"\t"pops"\t"af
#  }
#}' > /scratch/aob2x/dest/geo_endemic/geo_endemic.SNAPE.${maf}.${chr}.delim
#
