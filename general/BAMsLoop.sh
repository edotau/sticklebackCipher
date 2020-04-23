#!/bin/sh
set -e
for i in *.bam
do	
	sbatch $1 $i
done
exit 0
