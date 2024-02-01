#!/bin/bash
set -eo pipefail

REFERENCE=$1
COHORT_VCF=$2

# Ensure tools are available
export PATH=$PATH:$GATK_PATH

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <reference.fa> <cohort.genotypes.vcf.gz>"
    exit 1
fi

PREFIX=$(basename "$COHORT_VCF" | sed 's/\.vcf\.gz$//;s/\.vcf$//')

# Example for SNP part, similar logic can be applied for INDELs
RAW_SNP="${PREFIX}.raw.SNPs.vcf.gz"
FILTERED_SNP="${PREFIX}.filtered.SNPs.vcf.gz"

# Using named pipe for intermediate steps
mkfifo raw.snps.tmp default.snps.tmp frequency.snps.tmp

gatk --java-options "-Xmx4g" SelectVariants -R "$REFERENCE" -V "$COHORT_VCF" -select-type SNP -O raw.snps.tmp &
gatk --java-options "-Xmx4g" VariantFiltration -R "$REFERENCE" -V raw.snps.tmp --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "default_snp_filter" -O default.snps.tmp &
gatk --java-options "-Xmx4g" VariantFiltration -R "$REFERENCE" -V default.snps.tmp --filter-expression "AC < 4" --filter-name "frequency_filter" -O frequency.snps.tmp &
gatk --java-options "-Xmx4g" SelectVariants -R "$REFERENCE" -V frequency.snps.tmp --exclude-filtered -O "$FILTERED_SNP"

# INDEL processing variables
RAW_INDEL="${PREFIX}.raw.INDELs.vcf.gz"

FILTERED_INDEL="${PREFIX}.filtered.INDELs.vcf.gz"

# Create named pipes for intermediate INDEL steps
mkfifo raw.indels.tmp default.indels.tmp frequency.indels.tmp

# Process INDEL variants
gatk --java-options "-Xmx4g" SelectVariants -R "$REFERENCE" -V "$COHORT_VCF" -select-type INDEL -O raw.indels.tmp &
gatk --java-options "-Xmx4g" VariantFiltration -R "$REFERENCE" -V raw.indels.tmp --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "default_indel_filter" -O default.indels.tmp &
gatk --java-options "-Xmx4g" VariantFiltration -R "$REFERENCE" -V default.indels.tmp --filter-expression "AC < 4" --filter-name "frequency_filter" -O frequency.indels.tmp &
gatk --java-options "-Xmx4g" SelectVariants -R "$REFERENCE" -V frequency.indels.tmp --exclude-filtered -O "$FILTERED_INDEL"

# Cleanup named pipes
rm raw.snps.tmp default.snps.tmp frequency.snps.tmp
rm raw.indels.tmp default.indels.tmp frequency.indels.tmp
