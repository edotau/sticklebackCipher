#!/bin/sh -e
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --job-name=SNPs
module add GATK/4.1.0.0-gcb01
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)
rawSNP=${PREFIX}rawSNP.vcf.gz
defaultSNP=${PREFIX}.SNP.vcf.gz
filterSNP=${PREFIX}.Filtered.SNP.vcf.gz
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $VCF -select-type SNP -O $rawSNP
gatk --java-options "-Xmx16g" VariantFiltration -R $REF -V $rawSNP --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "default_snp_filter" --cluster-window-size 35 --cluster-size 3 -O $defaultSNP
gatk --java-options "-Xmx16g" SelectVariants -R $REF -V $defaultSNP --exclude-filtered -O $filterSNP
rm $defaultSNP ${defaultSNP}.tbi
