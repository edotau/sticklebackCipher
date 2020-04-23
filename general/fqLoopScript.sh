#!/bin/sh
set -e
for i in *R1*.{fastq,fq}*.gz
do	
	READ1=$i
	PREFIX=$(echo $READ1 | rev | cut -d '_' -f 2- | rev)
	SUFFIX=$(echo $READ1 | rev | cut -d '_' -f 1,2 | rev)
	READ2=$(echo $READ1 | rev | cut -d '_' -f 2- | rev).$SUFFIX
	
	sbatch $1 $READ1 $READ2
done
exit 0
