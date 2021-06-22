#!/bin/bash

for i in $(cat list_files.txt);
do
	for j in $(cat list_genomes.txt);
	do
		printf "%s " $j  >> ${i}/${i}.genomes_mapped_reads_coverage.txt; 
		tail ${i}/${i}_meandepth_original.bam.txt -n +2 | grep $j | awk '{print $4}' | paste -sd+ | bc >> ${i}/${i}.genomes_mapped_reads_coverage.txt;
	done
done
