#!/bin/bash

check_exit_status () {
  if [ ! "$2" -eq "0" ]; then
    echo "Step $1 failed with exit status $2"
    exit $2
  fi
  echo "Checked step $1"
}

bam="test_file"
sample="default_sample_name"
output="."
mean_read_length=75
min_BQ=25

# Credit: https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -mrl|--min-read-length)
    mean_read_length="$2"
    shift # past argument
    shift # past value
    ;;
    -mq|--min-BQ)
    min_BQ="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--help)
    echo "Usage:"
    echo "  weGot_thisCoverage [options] -h          Display this help message."
    echo "  weGot_thisCoverage [options] <bam_file_path> <sample_name> <output_dir>"
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

if [ $# != 3 ]
  then
    echo "ERROR: Need to supply <bam_file_path> <sample_name> <output_dir> as positional arguments, and no others"
    exit 1
fi

bam=$1; shift
sample=$1; shift
output=$1; shift

if [ ! -f "$bam" ]; then
  echo "ERROR: $bam does not exist"
  exit 1
fi

if [ ! -d $output/$sample/ ]; then
  mkdir -p $output/$sample/
fi

samtools coverage --min-read-len $mean_read_length --min-BQ $min_BQ -o $output/$sample/${sample}_meandepth.txt $bam

check_exit_status "samtools_coverage" $?

echo "BAM file: $bam" >> $output/$sample/${sample}.coverage.parameters.txt
echo "Sample name: $sample" >> $output/$sample/${sample}.coverage.parameters.txt
echo "Output directory: $output" >> $output/$sample/${sample}.coverage.parameters.txt
echo "Mean read length: $mean_read_length" >> $output/$sample/${sample}.coverage.parameters.txt
echo "Mean BQ: $min_BQ" >> $output/$sample/${sample}.coverage.parameters.txt

