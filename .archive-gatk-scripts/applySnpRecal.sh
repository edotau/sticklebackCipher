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
#indelDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDbLITCxMATAGenotypeGVCFs.Filtered.INDEL.vcf.gz
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)

#NOTE: second known site argument is to mask certain varience
gatk --java-options "-Xmx30g" VariantRecalibrator -R $REF -V $VCF \
	-resource:trainingSnpDb,known=false,training=true,truth=true,prior=8.0 $snpDb \
	-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode SNP \
	-O ${PREFIX}.recalibrated.SNPs.vcf.gz \
	--tranches-file ${PREFIX}.output.tranches --rscript-file ${PREFIX}.output.plots.R
