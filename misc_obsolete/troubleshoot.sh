#!/usr/bin/env bash
#
#SBATCH -J touchup # A single job name for the array
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --cpus-per-task=1
#SBATCH --mem 8G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/touchup.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/touchup.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/9.2.0  mvapich2/2.3.3 python/3.7.7
module load parallel

#python /scratch/aob2x/dest/DEST/mappingPipeline/scripts/PickleRef.py \
#  --ref /scratch/aob2x/dest/referenceGenome/r6/holo_dmel_6.12.fa \
#  --output /scratch/aob2x/dest/referenceGenome/r6/holo_dmel_6.12.fa.pickled.ref
#

### run as
sbatch --array=1-5 /scrach/aob2x/dest/DEST/misc_obsolete/troubleshoot.sh

### bad
#PA_li_09_fall
#WI_cp_14_fall
#TR_Yes_15_7

### good
#AT_Mau_14_02
#WI_cp_14_spring

#SLURM_ARRAY_TASK_ID=1
samps="PA_li_09_fall,WI_cp_14_fall,TR_Yes_15_7,AT_Mau_14_02,WI_cp_14_spring"

sample=$( echo ${samps} | tr ',' '\n' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

base_quality_threshold=15
illumina_quality_coding=1.8
minIndel=5
output=/project/berglandlab/DEST/dest_mapped/pipeline_output

zcat $output/${sample}/${sample}.mel_mpileup.txt.gz > /project/berglandlab/DEST/old/debug/${sample}.mel_mpileup.txt

python3 /scratch/aob2x/dest/DEST/mappingPipeline/scripts/Mpileup2Sync.py \
--mpileup /project/berglandlab/DEST/old/debug/${sample}.mel_mpileup.txt \
--ref /scratch/aob2x/dest/referenceGenome/r6/holo_dmel_6.12.fa.pickled.ref.ref \
--output /project/berglandlab/DEST/old/debug/${sample}.mel_mpileup.txt \
--base-quality-threshold $base_quality_threshold \
--coding $illumina_quality_coding \
--minIndel $minIndel
