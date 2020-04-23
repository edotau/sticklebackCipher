#!/bin/sh

gVcfMerged=stickleback.cohort

for i in *R1*.{fastq,fq}*
do
	PREFIX=$(basename $i _1.fq.gz)
	fqToBam=$(sbatch ./fqToBam.sh $i ${PREFIX}_2.fq.gz)
	dupsMates=$(sbatch --dependency=afterok:$fqToBam ./bamToGenotypeVcf.sh ${PREFIX}.bam)
	mergeGvcf=$(sbatch --dependency=singleton --job-name=$dupsMates ./combineGVCF.sh $gVcfMerged)
	genotypeVcf=$(sbatch --dependency=$mergeGvcf ./toGenotypeVcf ${gVcfMerged}.g.vcf.gz)
	sbatch ./getSNPs ${gVcfMerged}.vcf
	sbatch ./getIndels ${gVcfMerged}.vcf
done
exit 0
