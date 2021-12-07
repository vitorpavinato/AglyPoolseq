#!/usr/bin/env bash
#SBATCH --job-name=fastq_sanity_check # A single job name for the array
##SBATCH --ntasks-per-node=10 # one core
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00 ### most jobs should run in 60 minutes or less; the mitochondria takes a lot longer to run through pool-snp
#SBATCH --mem 20G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/fastqSanityCheck.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/fastqSanityCheck.%A_%a.err # Standard error
#SBATCH --account PAS1715

### modules
  python/3.6-conda5.2

### conda packages
source activate local_py3

### define a few things
  wd=/fs/scratch/PAS1715/aphidpool
  
### get job number
  # For the technical replicated run files
  prx=$( cat ${wd}/AglyPoolseq/populationInfo/fieldPools.csv | cut -f1,13 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )
  
  echo ${SLURM_ARRAY_TASK_ID}
  echo $pop
  echo $prx
  
  # For the replicated run files
  touch /fs/scratch/PAS1715/aphidpool/fastq/${prx}.R1.fastq.gz
  touch /fs/scratch/PAS1715/aphidpool/fastq/${prx}.R2.fastq.gz

### Run Fastq Sanity Check
${wd}/AglyPoolseq/mappingPipeline/scripts/fastq_sanity_check.sh \
/fs/scratch/PAS1715/aphidpool/fastq/${prx}.R1.fastq.gz

