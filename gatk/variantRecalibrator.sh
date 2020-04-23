#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=atacSNPs
module load GATK/4.1.3.0-gcb01
module load R gnuplot/4.6.4-fasrc02
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
MATAxLITC=/data/lowelab/edotau/sticklebackCipher/gatk/LITCxMATA.raw.truthSNPS.sitesonly.vcf.gz
atac=/data/lowelab/edotau/sticklebackCipher/gatk/LITCxMATA.atac.truthSNPs.sitesonly.vcf.gz
snpDb=/data/lowelab/edotau/sticklebackCipher/gatk/RABSxMatanuska.atacCohort.1stPass.SNPs.sitesonly.vcf.gz
#indelDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDbLITCxMATAGenotypeGVCFs.Filtered.INDEL.vcf.gz
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)


gatk --java-options "-Xmx16g" VariantRecalibrator -R $REF -V $VCF \
	-resource:atacHets,known=false,training=true,truth=true,prior=9.6 $atac \
	-resource:filteredMATAxLITC,known=false,training=true,truth=true,prior=9.5 $MATAxLITC \
	-resource:MxFw,known=false,training=true,truth=false,prior=8.0 $snpDb \
	-an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -mode SNP -AS \
	-O ${PREFIX}.VQSR.recalib.SNP.vcf.gz \
	--tranches-file ${PREFIX}.output.tranches --rscript-file ${PREFIX}.output.plots.R
