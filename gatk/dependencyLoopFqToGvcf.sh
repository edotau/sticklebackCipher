#!/bin/sh
set -e
gVcfMerged=stickleback.cohort

for i in *R1*.{fastq,fq}*.gz
do	
	READ1=$i
	SUFFIX=$(echo $READ1 | rev | cut -d '_' -f 1,2 | rev)
	READ2=$(echo $READ1 | rev | cut -d '_' -f 2- | rev).$SUFFIX
	
	fqToBam=$(sbatch ./fqToBam.sh $READ1 $READ2)
	dupsMates=$(sbatch --dependency=afterok:$fqToBam ./bamToGenotypeVcf.sh ${PREFIX}.bam)
done

mergeGvcf=$(sbatch --dependency=singleton --job-name=toGenotypeVcf ./combineGVCF.sh $gVcfMerged)
genotypeVcf=$(sbatch --dependency=$mergeGvcf ./toGenotypeVcf ${gVcfMerged}.g.vcf.gz)
sbatch ./getSNPs ${gVcfMerged}.vcf
sbatch ./getIndels ${gVcfMerged}.vcf
	
exit 0
