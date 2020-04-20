#!/bin/sh
for i in *_1.fq.gz
do
	PREFIX=$(basename $i _1.fq.gz)
	fqToBam=$(sbatch ./fqToBam.sh $i ${PREFIX}_2.fq.gz)
	dupsMates=$(sbatch --dependency=afterok:$fqToBam ./bamToGenotypeVcf.sh ${PREFIX}.bam)
done
exit 0
