#!/bin/sh

#  fastqc_sanity_check.sh
#  
#  Modified from https://sites.psu.edu/biomonika/2015/07/03/trim-now-or-never-pe-and-mp-adaptor-removal/
#  Created by Correa Pavinato, Vitor Antonio on 7/17/21.
#

function check_data {
R1=$1;
R2=$(echo $R1 | sed 's/R1.fastq.gz/R2.fastq.gz/g')

echo "Check if length of sequence equals quality score length"
bioawk -c fastx '{if (length($seq)!=length($qual)) print "Offending: " NR}' $R1
bioawk -c fastx '{if (length($seq)!=length($qual)) print "Offending: " NR}' $R2

echo "Finished validation."
echo "Check md5 sums:";
echo 'md5sum $R1'
echo 'md5sum $R2'
echo "——–"
}

#validate Fastq file
echo "File: " $1
check_data $1
