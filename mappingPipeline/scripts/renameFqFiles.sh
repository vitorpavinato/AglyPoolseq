#!/bin/bash

# usage: renameFqFile.sh <<list.files.prefix>>
# the list.files.prefix should only contain the file name withou .fastq.gz

for i in $(cat $1);
do 
	TMPNAME=`echo $i | awk '{split($0,a,"_"); print a[1]"_"a[3]"."a[2]}'`; 
	#echo $TMPNAME;
    #echo $i "\t" $TMPNAME
    mv ${i}.fastq.gz ${TMPNAME}.fastq.gz
done
