#!/usr/bin/env bash
#SBATCH -J run_baypass # A single job name for the array
#SBATCH -c 8
#SBATCH -N 1 # on one node
#SBATCH -t 168:00:00
#SBATCH --mem 4G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/run_baypass.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/run_baypass.%A_%a.err # Standard error
#SBATCH --account PAS1715

## Load modules
module load baypass

datafile="aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021"

wd="/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/baypass_poolsnp"

# Run BAYPASS CORE MODEL
# USED SEEDS 4040
#echo "running baypass core model"
#baypass -gfile ${wd}/${datafile}.genobaypass \
#        -poolsizefile ${wd}/${datafile}.poolsize -seed 4040 \
#        -outprefix ${wd}/core0 \
#        -nthreads 8 -nval 10000 -thin 25 -burnin 50000 -npilot 25 -pilotlength 1000;

# Run BAYPASS STD-IS MODEL
# USED SEEDS 6002, 3001, 5005
#echo "running baypass core std IS  model"
#baypass -gfile ${wd}/${datafile}.genobaypass \
#        -efile ${wd}/cov.biotypes -scalecov \
#        -poolsizefile ${wd}/${datafile}.poolsize -seed 3001 \
#        -outprefix ${wd}/stdis1 \
#        -nthreads 8 -nval 10000 -thin 25 -burnin 50000 -npilot 25 -pilotlength 1000;

# Run BAYPASS AUX MODEL
echo "running baypass AUX model"
baypass -gfile ${wd}/${datafile}.genobaypass \
        -efile ${wd}/cov.biotypes -scalecov -auxmodel \
        -omegafile ${wd}/core0/core0_mat_omega.out \
        -poolsizefile ${wd}/${datafile}.poolsize -seed 5001 \
        -outprefix ${wd}/aux0 \
        -nthreads 8 -nval 10000 -thin 25 -burnin 50000 -npilot 25 -pilotlength 500;


echo "done"
