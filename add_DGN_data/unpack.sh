#!/usr/bin/env bash
#
#SBATCH -J wide2long # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/wide2long.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/wide2long.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### unpack datasets
  #SLURM_ARRAY_TASK_ID=6

  name=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f2 -d' ' )
  url=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f3 -d' ' )
  fileName=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | rev | cut -f1 -d'/' | rev )

  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/wideData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

cd /scratch/aob2x/dest/dgn/wideData

if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then
  echo ${SLURM_ARRAY_TASK_ID}
  echo "foo1"
elif [ "${SLURM_ARRAY_TASK_ID}" == "2" ]; then
  echo "foo2"
elif [ "${SLURM_ARRAY_TASK_ID}" == "3" ]; then
  echo "foo3"
elif [ "${SLURM_ARRAY_TASK_ID}" == "5" ]; then
  echo "foo4"
elif [ "${SLURM_ARRAY_TASK_ID}" == "6" ]; then

  cd /scratch/aob2x/dest/dgn/wideData/NUZHDIN_sequences

  tar -xvf NUZHDIN_Chr2L_sequences.tar
  tar -xvf NUZHDIN_Chr2R_sequences.tar
  tar -xvf NUZHDIN_Chr3L_sequences.tar
  tar -xvf NUZHDIN_Chr3R_sequences.tar
  tar -xvf NUZHDIN_ChrX_sequences.tar


elif [ "${SLURM_ARRAY_TASK_ID}" == "7" ]; then
  echo "foo6"
elif [ "${SLURM_ARRAY_TASK_ID}" == "8" ]; then

  ###NONE!
  
fi
