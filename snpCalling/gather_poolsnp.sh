#!/usr/bin/env bash
#
#SBATCH -J gather_vcfs # A single job name for the array
##SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -c 2
#SBATCH -t 2:00:00 ### 6 hours
#SBATCH --mem 8G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/gather_vcfs.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/gather_vcfs.%A_%a.err # Standard error
#SBATCH --account PAS1715

### run as: sbatch --array=1-$( cat ${wd}/poolSNP_jobs.sample.csv | cut -f1 -d',' | sort | uniq | awk '{print NR}' | tail -n1 ) ${wd}/DEST/snpCalling/gather_poolsnp_paramtest.sh
### sacct -j 13029741
### cat /scratch/aob2x/dest/slurmOutput/split_and_run.12825614
module load htslib
module load bcftools/1.9.2
module load parallel2/19.10
 
wd="/fs/scratch/PAS1715/aphidpool"
outdir="/fs/scratch/PAS1715/aphidpool/sub_vcfs"

## Check if outdir exists
[ ! -d /fs/scratch/PAS1715/aphidpool/sub_bcf/ ] && mkdir /fs/scratch/PAS1715/aphidpool/sub_bcf/

#SLURM_ARRAY_TASK_ID=2
chr=$( cat ${wd}/poolSNP_jobs.csv | cut -f1 -d',' | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d" )

popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
#maf=01; mac=50; popSet="all"; method="PoolSNP"; version="paramTest"

ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz | sort -t"_" -k2,2 -k3g,3  | \
#grep /${chr}_ > /fs/scratch/PAS1715/aphidpool/sub_vcfs/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort
grep /${chr}_

#bcftools concat \
#--threads 20 \
#-f /fs/scratch/PAS1715/aphidpool/sub_vcfs/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
#-O v \
#-n \
#-o /fs/scratch/PAS1715/aphidpool/sub_bcf/aphidpool.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.bcf
