#!/bin/bash
set -eo pipefail

# Ensure tools are available
export PATH=$PATH:$GATK_PATH

# Assign arguments to variables
BAM=$1
REFERENCE=$2
SNP_VCF=$3
INDEL_VCF=$4

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $(basename "$0") <reference.fa> <sample.bam> <snp.vcf> <indel.vcf>"
    exit 1
fi

# Directories and file names
DIR="bqsr-genotype-vcfs"
RECALIBRATED_BAMS_DIR="recalibrated-bams"
mkdir -p "$DIR" "$RECALIBRATED_BAMS_DIR"

PREFIX=$(basename "$BAM" .bam)
BAM_RECALIBRATED="$RECALIBRATED_BAMS_DIR/${PREFIX}.recal.bam"
RECALIBRATED_TABLE="$RECALIBRATED_BAMS_DIR/${PREFIX}.recal_data.table"
GENOTYPE_VCF="$DIR/${PREFIX}.g.vcf.gz"

# Step 1: Base quality score recalibration
gatk BaseRecalibrator -R "$REFERENCE" \
    --known-sites "$SNP_VCF" \
    --known-sites "$INDEL_VCF" \
    -I "$BAM" \
    -O "$RECALIBRATED_TABLE"

# Step 2: Apply BQSR
gatk ApplyBQSR -R "$REFERENCE" \
    -I "$BAM" \
    --bqsr-recal-file "$RECALIBRATED_TABLE" \
    -O "$BAM_RECALIBRATED"

# Step 3: Index recalibrated BAM
samtools index "$BAM_RECALIBRATED"

# Step 4: Variant calling with HaplotypeCaller
gatk HaplotypeCaller \
    --input "$BAM_RECALIBRATED" \
    --output "$GENOTYPE_VCF" \
    --reference "$REFERENCE" \
    -ERC GVCF
