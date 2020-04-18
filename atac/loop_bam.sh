#!/bin/sh
for i in *.sorted
do
	sbatch --mem=32G --nodes=1 --ntasks=1 --cpus-per-task=8 $1 $i
done
exit 0
