# Scripts to download, map, call polymorphism in pooled sequencing data-sets for Drosophila

## Description
> This set of scripts provides a pipeline to build wholeGenomeSync files for each population sample from raw FASTQ data

### 0. Define working directory
```bash
   wd=/scratch/aob2x/dest
```

### 1. Download data from SRA
```bash
   sbatch --array=1-$( wc -l < ${wd}/DEST/populationInfo/samps.csv ) \
   ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
```
A few timed out, restart. <br/>
``` bash
  sbatch --array=$( sacct -u aob2x -j 9199377 | grep "TIMEOUT" | cut -f1 -d' ' | cut -f2 -d'_' | tr '\n' ',' ) \
  ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
```
sacct -u aob2x -j 9213810

A few more timed out again, restart those. Now up to 36 hour time limit. Plus Maine samples get fixed.
``` bash
  sbatch --array=$( sacct -u aob2x -j 9213810 | grep "TIMEOUT" | cut -f1 -d' ' | cut -f2 -d'_' | tr '\n' ',' ) \
  ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
```
sacct -u aob2x -j 9248781


### 2. Check that data are in the correct FASTQ format
Double check that all downloaded data are in Fastq33. Uses script from [here](https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py). </br>
```bash
  sbatch ${wd}/DEST/mappingPipeline/scripts/check_fastq_encoding.sh

  grep -v "1.8" ${wd}/fastq/qualEncodings.delim
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
