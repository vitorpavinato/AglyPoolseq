#!/usr/bin/env bash
#SBATCH -J run_selestim_pods # A single job name for the array
#SBATCH -c 28 # on one entire node
#SBATCH -N 1  # on one node
#SBATCH -t 168:00:00
#SBATCH --mem 4G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/run_selestim_pods.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/run_selestim_pods.%A_%a.err # Standard error
#SBATCH --account PAS1715

## Load modules
module load SelEstim

datafile="aphidpool.PoolSeq.PoolSNP.05.5.24Jun2021_rand"

wd="/fs/scratch/PAS1715/aphidpool/results/aggregated_data/minmaxcov_4_99/selestim_poolsnp"

# Run SelEstim
echo "running selestim"
selestim -file ${wd}/${datafile}.genoselestim \
         -npilot 20 -lpilot 500 -burnin  120000 -length 50000 -thin 50 -seed 1234 -pool \
         -outputs ${wd}/output_2/ \
         -calibration_only

echo "done"




