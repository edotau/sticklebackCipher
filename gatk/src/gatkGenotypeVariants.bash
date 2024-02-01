#!/bin/bash
set -eo pipefail

# Check the number of arguments
if [ "$#" -ne 2 ]; then
    echo "
    Usage:
        $(basename "$0") <reference.fa> <raw.genotypes.vcf.gz>
    "
    exit 1
fi

# Assign arguments to variables
REFERENCE=$1
COHORT_VCF=$2

# Ensure tools are available
export PATH=$PATH:$GATK_PATH

# Define file paths and prefixes
PREFIX=$(basename "$COHORT_VCF" | sed 's/\.vcf\.gz$//;s/\.vcf$//')
SNPS_DB=/path/to/trainingSnpDb.vcf.gz  # Update this path to your SNP database
RAW_SNP="${PREFIX}.raw.SNPs.vcf.gz"

BQSR_SNP="${PREFIX}.bqsr"
TRACE_FILE="${BQSR_SNP}.output.tranches"

# Select SNPs from the cohort VCF
gatk --java-options "-Xmx16g" SelectVariants \
    -R "$REFERENCE" \
    -V "$COHORT_VCF" \
    -select-type SNP \
    -O "$RAW_SNP"

# Recalibrate variant quality scores for SNPs
gatk --java-options "-Xmx30g" VariantRecalibrator \
    -R "$REFERENCE" \
    -V "$RAW_SNP" \
    -resource:trainingSnpDb,known=false,training=true,truth=true,prior=8.0 "$SNPS_DB" \
    -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \
    -mode SNP \
    -O "$BQSR_SNP" \
    --tranches-file "$TRACE_FILE" \
    --rscript-file "${BQSR_SNP}.output.plots.R"

# Apply the recalibration to the SNP set
gatk --java-options "-Xmx30g" ApplyVQSR \
    -R "$REFERENCE" \
    -V "$RAW_SNP" \
    -ts-filter-level 99.9 \
    -recal-file "$BQSR_SNP" \
    --tranches-file "$TRACE_FILE" \
    -O "${PREFIX}.vqsr.SNPs.vcf.gz"

# Filter the final set of SNPs
gatk --java-options "-Xmx8g" VariantFiltration \
    -V "${PREFIX}.vqsr.SNPs.vcf.gz" \
    -O "${PREFIX}.vqsr.ase.genotype.SNPs.vcf.gz" \
    --cluster-window-size 35 \
    --cluster-size 3 \
    --filter-name "FS" --filter-expression "FS > 30.0" \
    --filter-name "QD" --filter-expression "QD < 2.0"
