#!/usr/bin/env bash
#
#SBATCH -J reheader # A single job name for the array
#SBATCH --ntasks-per-node=9 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/reheader.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/reheader.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/dest/DEST/utils/fix_Ukraine_sampleNames/remakeGDS.sh
### sacct -j 19747891
### cat /scratch/aob2x/dest/slurmOutput/reheader.19747891

module load intel/18.0 intelmpi/18.0 R/3.6.3 bcftools parallel

wd="/project/berglandlab/DEST"

### VCFs
  vcf_file_1="vcf/dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
  vcf_file_2="vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
  vcf_file_3="vcf/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz"

  makeGDS () {
    f1=${1} # f1=vcf/dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz
    f2=$( echo ${f1} | sed 's/\.gz//g' | sed 's/^vcf/gds/g' )
    wd="/project/berglandlab/DEST"


    bcftools view \
    -O v \
    --threads 3 \
    ${wd}/${f1} > \
    /scratch/aob2x/${f2}

    Rscript /scratch/aob2x/dest/DEST/snpCalling/vcf2gds.R ${f2}

  }
  export -f makeGDS

  parallel --jobs 3 makeGDS ::: $vcf_file_1 $vcf_file_2 $vcf_file_3
