#!/bin/bash
## --ntasks-per-node can be increased upto 48 on Pitzer
#SBATCH --nodes=1 --ntasks-per-node=28 
#SBATCH --time=05:00:00
#SBATCH --job-name Blast
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/run_blast.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/run_blast.%A_%a.err # Standard error
#SBATCH --account=PAS1715

module load blast/2.11.0+
module load blast-database/2021-05

## Set the running environment
# PATHS
wd=/fs/scratch/PAS1715/aphidpool
inputDir=${wd}/results/aggregated_data/minmaxcov_4_99/signf_genes_annotation
outputDir=${inputDir}

# COPY FILES TO TMPDIR
cp ${inputDir}/outliers_genes_bedtools-snpeff-annR-windows_aa.fasta $TMPDIR
cd $TMPDIR

## Run blastp
echo "running blastp"
blastp -db nr -query outliers_genes_bedtools-snpeff-annR-windows_aa.fasta -num_threads 16 -evalue 1e-6 -outfmt 5 -out outliers_genes_bedtools-snpeff-annR-windows_blastp.out

## Move the results back
cp outliers_genes_bedtools-snpeff-annR-windows_blastp.out $SLURM_SUBMIT_DIR/${outputDir}

echo "done"
