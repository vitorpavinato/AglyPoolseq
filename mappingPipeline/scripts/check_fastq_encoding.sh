#!/bin/bash
#SBATCH --job-name=check_fastq_encoding
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem 1G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/check_fastq_encoding.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/check_fastq_encoding.%A_%a.err # Standard error
#SBATCH --account PAS1715

wd=/fs/scratch/PAS1715/aphidpool

chmod +x ${wd}/AglyPoolseq/mappingPipeline/scripts/guess_encoding.py

if test -f ${wd}/fastq/qualEncodings.delim; then
    rm ${wd}/fastq/qualEncodings.delim
    touch ${wd}/fastq/qualEncodings.delim
else
    touch ${wd}/fastq/qualEncodings.delim
fi

for f in ${wd}/fastq/*fastq.gz; do
  # f=${wd}/fastq/Pool7_AddRun.R1.fastq.gz
  zcat $f | head -n 5000 | awk 'NR % 4 == 0' | \
  ${wd}/AglyPoolseq/mappingPipeline/scripts/guess_encoding.py | grep -v "#" | sed -E "s|^|${f}\t|g" >> ${wd}/fastq/qualEncodings.delim
done

sleep 300
