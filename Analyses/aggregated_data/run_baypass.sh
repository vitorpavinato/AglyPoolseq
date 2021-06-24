#!/usr/bin/env bash
#SBATCH -J run_baypass # A single job name for the array
##SBATCH --ntasks-per-node=1 # one core
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 20:00:00 ### most jobs should run in 60 minutes or less; the mitochondria takes a lot longer to run through pool-snp
#SBATCH --mem 10G
#SBATCH -o /fs/scratch/PAS1715/aphidpool/slurmOutput/run_baypass.%A_%a.out # Standard output
#SBATCH -e /fs/scratch/PAS1715/aphidpool/slurmOutput/run_baypass.%A_%a.err # Standard error
#SBATCH --account PAS1715

## Load modules
module load baypass

popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
missing=${6}
pools=${7}
#maf=001; mac=100; popSet="all"; method="PoolSNP"; version="paramTest"; missing=05; pools="allpools"

wd="/fs/scratch/PAS1715/aphidpool"

# Run BAYPASS
echo "running baypass"
baypass -gfile ${wd}/data/aggregated_data/genobaypass/aphidpool.${popSet}.${method}.${maf}.${mac}.${version}.vcf.missfrac${missing}.${pools}.genobaypass \
        -poolsizefile ${wd}/data/aggregated_data/genobaypass/aphidpool.${popSet}.${method}.${maf}.${mac}.${version}.vcf.missfrac${missing}.${pools}.poolsize \
        -outprefix core.${popSet}.${method}.${maf}.${mac}.${version}.missfrac${missing}.${pools} -nthreads 5 -pilotlength 500 -burnin 5000 -seed 5001;

echo "done"
