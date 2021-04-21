#!/bin/bash

check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2"
    exit $2
  fi
  echo "Checked step $1"
}

suffix="test_file"
input="."
sample="default_sample_name"
output="."

# Credit: https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo "Usage:"
    echo "  aggregate_bams -h Display this help message."
    echo "  aggregate_bams <bam_suffix_file> <input_dir> <sample_name> <output_dir>"
    exit 0
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") #save it to an array
    shift
    ;;
esac
done

set -- "${POSITIONAL[@]}"

if [ $# != 4 ]
  then
    echo "ERROR: Need to supply <bam_suffix_file> <input_dir> <sample_name> <output_dir> as positional arguments, and no others"
    exit 1
fi

suffix=$1; shift
input=$1; shift
sample=$1; shift
output=$1; shift

if [ ! -d $output/$sample/ ]; then
  mkdir -p $output/$sample/
fi

# create a list of files
ls -1 $input/${sample}_*/${sample}_*.${suffix}.bam > $output/$sample/${sample}_bamlist.txt

BAMLIST=$output/$sample/${sample}_bamlist.txt
NUM_FILES=$(wc -l < "$BAMLIST")

if [ "${NUM_FILES}" == "5" ]; then
    echo "There are 5"
    java -jar $PICARD MergeSamFiles I=$input/${sample}_140711/${sample}_140711.${suffix}.bam I=$input/${sample}_140722/${sample}_140722.${suffix}.bam I=$input/${sample}_140725/${sample}_140725.${suffix}.bam I=$input/${sample}_140730/${sample}_140730.${suffix}.bam I=$input/${sample}_AddRun/${sample}_AddRun.${suffix}.bam SO=coordinate USE_THREADING=true O=$output/$sample/${sample}.aggregated.bam

else
    echo "There are 4"
    java -jar $PICARD MergeSamFiles I=$input/${sample}_140711/${sample}_140711.${suffix}.bam I=$input/${sample}_140722/${sample}_140722.${suffix}.bam I=$input/${sample}_140725/${sample}_140725.${suffix}.bam I=$input/${sample}_140730/${sample}_140730.${suffix}.bam SO=coordinate USE_THREADING=true O=$output/$sample/${sample}.aggregated.bam

fi
    
check_exit_status "Picard_MergeSamFiles" $?
  
  
