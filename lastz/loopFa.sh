#!/bin/sh -e
DIR=$1
SCRIPT=$2
QUERY=$3
for i in ${DIR}/*.fa
do
        sbatch --mem=8G --ntasks=1 --cpus-per-task=4 $SCRIPT $i $QUERY
done
