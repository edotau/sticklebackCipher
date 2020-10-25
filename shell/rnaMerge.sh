#!/bin/sh

for i in *_01_dupMarked.bam
do
	SAMPLE=$(basename $i _01_dupMarked.bam)
	sbatch ./run_pipeline.sh $i ${SAMPLE}_02_dupMarked.bam ${SAMPLE}.bam
done
exit 0
