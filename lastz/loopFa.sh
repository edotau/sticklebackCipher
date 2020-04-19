#!/bin/sh -e
DIR=${1%/}
SCRIPT=$2
QUERY=$3
for i in ${DIR}/*.fa
do
        sbatch $SCRIPT $i $QUERY
done
