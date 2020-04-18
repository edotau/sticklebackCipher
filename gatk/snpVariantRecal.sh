#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --job-name=variantRecalibrator
module load GATK/4.1.3.0-gcb01
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
snpDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDb/LITCxMATAGenotypeGVCFsrawSNP.vcf.gz
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)

rawSNP=${PREFIX}.raw.SNPs.vcf.gz
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $VCF \
	-select-type SNP -O $rawSNP
#output files for variance recal
snpRecal=${PREFIX}.recal
tsFile=${PREFIX}.output.tranches
#NOTE: second known site argument is to mask certain variance
gatk --java-options "-Xmx30g" VariantRecalibrator -R $REF -V $rawSNP \
	-resource:trainingSnpDb,known=false,training=true,truth=true,prior=8.0 $snpDb \
	-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP \
	-O ${PREFIX}.recal --tranches-file $tsFile --rscript-file ${PREFIX}.output.plots.R
# Command 20: ApplyVQSR
#NOTE: second known site argument is to mask certain varience
gatk --java-options "-Xmx30g" ApplyVQSR -R $REF -V $rawSNP \
	-ts-filter-level 99.9 -recal-file $snpRecal --tranches-file $tsFile
	-O ${PREFIX}.recal.SNPs.vcf.gz

gatk --java-options "-Xmx8g" VariantFiltration -V ${PREFIX}.recal.SNPs.vcf.gz -O ${PREFIX}.ASE.SNP,final.vcf.gz \
--cluster-window-size 35 --cluster-size 3 --filter-name FS --filter-expression "FS > 30.0" \
--filter-name QD --filter-expression "QD < 2.0"