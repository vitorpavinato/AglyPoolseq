#!/usr/bin/env bash
#
#SBATCH -J unpack # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/unpack.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/unpack.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### unpack datasets
  #SLURM_ARRAY_TASK_ID=6

  name=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f2 -d' ' )
  url=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f3 -d' ' )
  fileName=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | rev | cut -f1 -d'/' | rev )


cd /scratch/aob2x/dest/dgn/wideData

if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then
  #tar -zvx \
  #--directory=/scratch/aob2x/dest/dgn/rawData \
  #--file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/
  tar -xvf dpgp2_Chr2L_sequences.tar
  tar -xvf dpgp2_Chr2R_sequences.tar
  tar -xvf dpgp2_Chr3L_sequences.tar
  tar -xvf dpgp2_Chr3R_sequences.tar
  tar -xvf dpgp2_ChrX_sequences.tar

  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

  #ls *.seq | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/DPGP2_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "2" ]; then
  tar -jvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

elif [ "${SLURM_ARRAY_TASK_ID}" == "3" ]; then
  #tar -jvx \
  #--directory=/scratch/aob2x/dest/dgn/rawData \
  #--file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/
  tar -xvf dgrp_Chr2L.tar
  tar -xvf dgrp_Chr2R.tar
  tar -xvf dgrp_Chr3L.tar
  tar -xvf dgrp_Chr3R.tar
  tar -xvf dgrp_ChrX.tar

  ls RAL* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/DGRP_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "5" ]; then
    tar -zvx \
    --directory=/scratch/aob2x/dest/dgn/rawData \
    --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences

    tar -xvf CLARK_Chr2L_sequences.tar
    cd /scratch/aob2x/dest/dgn/rawData/CLARK/CLARK_Chr2L
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/PEK_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/INY_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/NL_{}

    tar -xvf CLARK_Chr2R_sequences.tar
    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr2R
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/PEK_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/INY_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/NL_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/TAS_{}

    tar -xvf CLARK_Chr3L_sequences.tar
    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr3L
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/PEK_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/INY_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/NL_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/TAS_{}

    tar -xvf CLARK_Chr3R_sequences.tar
    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr3R
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/PEK_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/INY_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/NL_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/TAS_{}

    tar -xvf CLARK_ChrX_sequences.tar
    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_ChrX
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/PEK_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/INY_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/NL_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/TAS_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "6" ]; then
  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences

  tar -xvf NUZHDIN_Chr2L_sequences.tar
  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr2L
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/WCA_{}

  tar -xvf NUZHDIN_Chr2R_sequences.tar
  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr2R
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/WCA_{}

  tar -xvf NUZHDIN_Chr3L_sequences.tar
  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr3L
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/WCA_{}

  tar -xvf NUZHDIN_Chr3R_sequences.tar
  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr3R
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/WCA_{}

  tar -xvf NUZHDIN_ChrX_sequences.tar
  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_ChrX
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/WCA_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "7" ]; then
  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

elif [ "${SLURM_ARRAY_TASK_ID}" == "8" ]; then
  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/
  ls Simulans* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/SIM_{}
fi
