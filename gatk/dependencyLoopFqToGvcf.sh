#!/bin/sh
set -e
gVcfMerged=stickleback.cohort

BAMS="basic_alignment"
QCedBAMS="gatkValidBAMs"
toGenotypeVcf="toGenotypeVcfs"

mkdir -p $BAMS
mkdir -p $QCedBAMS
mkdir -p $toGenotypeVcf

bamToQC="bamToQC"
for i in *R1*.{fastq,fq}*.gz
do	
	READ1=$i
	PREFIX=$(echo $READ1 | rev | cut -d '_' -f 2- | rev)
	SUFFIX=$(echo $READ1 | rev | cut -d '_' -f 1,2 | rev)
	READ2=$(echo $READ1 | rev | cut -d '_' -f 2- | rev).$SUFFIX
	
	fqToBam=$(sbatch ./fqToBam.sh $READ1 $READ2 $BAMS)
	dupsMates=$(sbatch --dependency=afterok:$fqToBam ./basicBamToMarkdupsFixmate.sh $BAMS/${PREFIX}.bam $QCedBAMS)
	genotypeBAMs=$(sbatch --dependency=afterok:$dupsMates --job-name=$bamToQC gVcfGatk.sh $QCedBAMS/${PREFIX}.bam $toGenotypeVcf)
done

mergeGvcf=$(sbatch --dependency=singleton  --job-name=$bamToQC ./combineGVCF.sh $gVcfMerged $toGenotypeVcf)
genotypeVcf=$(sbatch --dependency=$mergeGvcf ./toGenotypeVcf ${gVcfMerged}.g.vcf.gz)
sbatch ./getSNPs ${gVcfMerged}.vcf
sbatch ./getIndels ${gVcfMerged}.vcf
	
exit 0
