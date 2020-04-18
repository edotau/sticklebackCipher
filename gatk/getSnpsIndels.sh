#!/bin/sh -e
module add GATK/4.1.3.0-gcb01
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)

rawSNP=${PREFIX}rawSNP.vcf.gz
defaultSNP=${PREFIX}.SNP.vcf.gz
frequencySNP=${PREFIX}.better.SNP.vcf.gz
filterSNP=${PREFIX}.Filtered.SNP.vcf.gz
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $VCF -select-type SNP -O $rawSNP
gatk --java-options "-Xmx16g" VariantFiltration -R $REF -V $rawSNP --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "default_snp_filter" -O $defaultSNP
gatk --java-options "-Xmx16g" VariantFiltration -R $REF -V $defaultSNP --filter-expression "AC < 4" --filter-name "frequency_filter" -O $frequencySNP
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $frequencySNP --exclude-filtered -O $filterSNP
rm $rawSNP $defaultSNP $frequencySNP ${rawSNP}.tbi ${defaultSNP}.tbi ${frequencySNP}.tbi

rawINDEL=${PREFIX}rawINDEL.vcf.gz
defaultINDEL=${PREFIX}.INDEL.vcf.gz
frequencyINDEL=${PREFIX}.better.INDEL.vcf.gz
filterINDEL=${PREFIX}.Filtered.INDEL.vcf.gz
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $VCF -select-type INDEL -O $rawINDEL
gatk --java-options "-Xmx16g" VariantFiltration -R $REF -V $rawINDEL --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "default_indel_filter" -O $defaultINDEL
gatk --java-options "-Xmx16g" VariantFiltration -R $REF -V $defaultINDEL --filter-expression "AC < 4" --filter-name "frequency_filter" -O $frequencyINDEL
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $frequencyINDEL -exclude-filtered -O $filterINDEL
rm $rawINDEL $defaultINDEL $frequencyINDEL ${rawINDEL}.tbi ${defaultINDEL}.tbi ${frequencyINDEL}.tbi

##Base Recalibration
#gatk --java-options "-Xmx16g" BaseRecalibrator -R $REF \
#	--known-sites \
#	--known-sites
