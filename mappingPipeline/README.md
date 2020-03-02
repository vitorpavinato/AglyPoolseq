# Scripts to download, map, call polymorphism in pooled sequencing data-sets for Drosophila

## Description
> This set of scripts provides a pipeline to build wholeGenomeSync files for each population sample from raw FASTQ data

### 0. Define working directory
```bash
   wd=/scratch/aob2x/dest
```

### 1. Download data from SRA
```bash
   sbatch --array=1-$( wc -l < ${wd}/DEST/populationInfo/samps.csv ) ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
```



















### notes for using Docker

  docker pull alanbergland/dest_mapping ### pull from dockerhub
  docker images  ### see what I have
  docker ps ### gives list of docker images running
  docker run -it NAME ### interactive
  docker run -v /PATH_TO_FASTQ/:/data/ NAME ./mapping.sh /data/fastq1 /data/fastq2


### for singularity on Rivanna
  module load singularity

  ### pull and convert to singularity file
  singularity pull -n dest_mapping.simg \
  docker://alanbergland/dest_mapping

  singularity exec dest_mapping.simg samtools -t ${SLURM_CPUS_PER_TASK}

  ### pull and run basic command
  singularity exec \
  docker://alanbergland/dest_mapping ls /opt



  singularity exec -c -B /PATH_TO_FASTQ:/data \
  docker://alanbergland/dest_mapping \
  /DIR/mapping.sh /data/fastq1.fq /data/fastq2.fq

  $( echo ${SLURM_CPUS_PER_TASK} )

  singularity exec -c -B /PATH_TO_FASTQ:/data \
  docker://alanbergland/dest_mapping \
  /DIR/mapping.sh /data/fastq1.fq /data/fastq2.fq ${SLURM_CPUS_PER_TASK}
