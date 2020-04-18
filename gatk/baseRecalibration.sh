#!/bin/sh
module load GATK/4.1.3.0-gcb01
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
snpDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDbLITCxMATAGenotypeGVCFs.Filtered.SNP.vcf.gz
indelDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDbLITCxMATAGenotypeGVCFs.Filtered.INDEL.vcf.gz
BAM=$1
PREFIX=$(basename $BAM .bam)

#NOTE: second known site argument is to mask certain varience

gatk --java-options "-Xmx16g" BaseRecalibrator -R $REF \
	--known-sites $snpDb \
       	--known-sites  $indelDb \
	-I $BAM -O ${PREFIX}.recal_data.table

gatk --java-options "-Xmx16g" ApplyBQSR -R $REF \
	-I $BAM --bqsr-recal-file ${PREFIX}.recal_data.table \
	-O ${PREFIX}.recal.bam

gatk HaplotypeCaller --java-options "-Xmx20G" \
	--input $BAM --output ${PREFIX}.wgs.g.vcf \
	--reference $REF \
	-ERC GVCF --native-pair-hmm-threads 8
