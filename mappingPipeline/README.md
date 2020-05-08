# Scripts to download, map, call polymorphism in pooled sequencing data-sets for Drosophila

## Description
> This set of scripts provides a pipeline to build wholeGenomeSync files for each population sample from raw FASTQ data and defines a Dockerfile to build a docker image which can act as a standalone tool to run the pipeline.

### 0. Define working directory
```bash
   wd=/scratch/aob2x/dest
```

### 1. Download data from SRA (specify 72 hour time limit)
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

One more!
``` bash
sbatch --array=175 ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
```
sacct -u aob2x -j 9248827

Ack! They keep timing out! now up to 72 hours.
``` bash
  sbatch --array=$( sacct -u aob2x -j 9248781 | grep "TIMEOUT" | cut -f1 -d' ' | cut -f2 -d'_' | tr '\n' ',' ) \
  ${wd}/DEST/mappingPipeline/scripts/downloadSRA.sh
```
sacct -u aob2x -j 9362168

### 2. Check that data are in the correct FASTQ format
Double check that all downloaded data are in Fastq33. Uses script from [here](https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py). </br>
```bash
  sbatch ${wd}/DEST/mappingPipeline/scripts/check_fastq_encoding.sh

  grep -v "1.8" ${wd}/fastq/qualEncodings.delim
```
### 3. Build singularity container from docker image
```bash
    module load singularity
    singularity pull docker://jho5ze/dmelsync:latest
```

### 4. Run the singularity container
```bash
    singularity run dmelsync_latest.sif <read_1> <read_2> <sample_name> <output_folder> <num_cores_to_use>
```
#### Input
* read_1 full path
* read_2 full path
* name of sample
* full path of output directory

#### Output contained in the output directory under a folder named for the sample name
* directory of fastq analysis for trimmed and untrimmed reads
* intervals file used for GATK IndelRealigner
* original bam file containing all reads (original.bam)
* simulans contaminants bam file (sim.bam)
* melanogaster reads (mel.bam)
* genomewide SYNC file (genomewide.sync)

### How to edit the pipeline

* ToDo: Create team on Docker to allow collaborative editing of containers

* Comment out the lines in the Dockerfile defining the ```ENTRYPOINT``` and ```CMD```

* Build the docker on your personal computer using ```docker build -t <image_name> .``` in the directory containing the Dockerfile (same as this README)

* Modify the pipeline script

* Run the docker using volumes (-v flag) to bind any desired directories (including the one in which you are modifying the pipeline script) to the image instance and -it to make it interactive: ```docker run -v <local_path>:<docker_path> -it <image_name>```

* Run the pipeline script, or any modifications you mean to make to the script, in this interactive version of the docker to debug, repeating the last two steps (or simply modifying the script within the interactive shell with vim) until you are satisfied with the changes.

* Push your changes to https://github.com/alanbergland/DEST.git

* Reassign the ARG CACHE_BUST line in the Dockerfile to the current month, day, year, and time so that the build process knows to rerun the following lines and repull the repo with the updated scripts.

* Uncomment the ```ENTRYPOINT``` and ```CMD``` lines

* Run the build process again, and then push to docker hub (currently you will have to make your own repo for this)

* Use ```singularity pull docker://<user>/<image_name>:<version>``` to get the updated image on your cluster
