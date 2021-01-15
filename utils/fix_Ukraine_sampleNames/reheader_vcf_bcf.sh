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

### sbatch /scratch/aob2x/dest/DEST/utils/fix_Ukraine_sampleNames/reheader_vcf_bcf.sh
### sacct -j 19602112
### cat /scratch/aob2x/dest/slurmOutput/reheader.19602112

module load bcftools parallel samtools

wd="/project/berglandlab/DEST"

### VCFs
  vcf_file_1="vcf/dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
  vcf_file_2="vcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz"
  vcf_file_3="vcf/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz"

  renameVCF () {
    f1=${1}

    echo ${f1}

    bcftools view -h ${wd}/${f1} |
    sed 's/UA_Pir_14_26/UA_Pyr_14_26/g' |
    sed 's/UA_Pir_15_21/UA_Pyr_15_21/g' |
    sed 's/UA_Pyr_16_48/UA_Pir_16_48/g' |
    sed 's/FLoat/Float/g'> \
    ${wd}/${f1}.header

    bcftools reheader \
    -h ${wd}/${f1}.header \
    ${wd}/${f1} > \
    ${wd}/${f1}.reheader

    mv ${wd}/${f1} ${wd}/old/${f1}.oldheader

    mv ${wd}/${f1}.reheader ${wd}/${f1}
    tabix -p vcf -f ${wd}/${f1}

    rm ${wd}/${f1}.header
  }
  export -f renameVCF

  parallel --jobs 3 renameVCF ::: $vcf_file_1 $vcf_file_2 $vcf_file_3

### BCFs
  bcf_file_1="bcf/dest.all.PoolSNP.001.50.10Nov2020.header.bcf"
  bcf_file_2="bcf/dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf"
  bcf_file_3="bcf/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.header.bcf"

  renameBCF () {
    f1=${1}
    # f1=$bcf_file_1
    echo ${f1}

    bcftools view -h ${wd}/${f1} |
    sed 's/UA_Pir_14_26/UA_Pyr_14_26/g' |
    sed 's/UA_Pir_15_21/UA_Pyr_15_21/g' |
    sed 's/UA_Pyr_16_48/UA_Pir_16_48/g' |
    sed 's/FLoat/Float/g'> \
    ${wd}/${f1}.header

    bcftools reheader \
    -h ${wd}/${f1}.header \
    ${wd}/${f1} |
    bcftools view \
    -O b --threads 3 > \
    ${wd}/${f1}.reheader

    mv ${wd}/${f1} ${wd}/old/${f1}.oldheader

    mv ${wd}/${f1}.reheader ${wd}/${f1}

    rm ${wd}/${f1}.header
  }
  export -f renameBCF

  parallel --jobs 3 renameBCF ::: $vcf_file_1 $vcf_file_2 $vcf_file_3
