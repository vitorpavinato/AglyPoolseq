#!/bin/sh

#SBATCH --job-name=docker_pull_osc
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --mem 4G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/docker_pull_osc.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/docker_pull_osc.%A_%a.err # Standard error
#SBATCH --account PAS1715

module load singularity

wd=/fs/scratch/PAS1715/aphidpool

cd ${wd}

singularity --debug pull docker://vpavinato/aglypoolseq:latest
