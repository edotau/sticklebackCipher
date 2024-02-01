#!/bin/bash
set -eo pipefail

REFERENCE=$1
COHORT_VCF=$2

# Ensure tools are available
export PATH=$PATH:$GATK_PATH

if [ "$#" -ne 2 ]; then
    echo "
        Usage: $(basename "$0") <reference.fa> <cohort.genotypes.vcf.gz>
    "
    exit 1
fi

# Validate required tools
for cmd in gatk; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in the current PATH."
        exit 1
    fi
done

PREFIX=$(basename "$COHORT_VCF" | sed 's/\.vcf\.gz$//;s/\.vcf$//')

# Example for SNP part, similar logic can be applied for INDELs
RAW_SNP="${PREFIX}.raw.SNPs.vcf.gz"
DEFAULT_SNP="${PREFIX}.default.SNPs.vcf.gz"
FREQUENCY_SNP="${PREFIX}.frequency.SNPs.vcf.gz"
FILTERED_SNP="${PREFIX}.filtered.SNPs.vcf.gz"

# Using named pipe for intermediate steps
mkfifo "$RAW_SNP" "$DEFAULT_SNP" "$FREQUENCY_SNP"

OS=$(uname -s)
case "$OS" in
    Linux*)
        MEM=$(awk '/MemAvailable/ {printf "%d", $2/1024/1024*0.9 }' /proc/meminfo)
        ;;
    Darwin*)
        MEM=$(sysctl -n hw.memsize | awk '{printf "%d", $1/1024/1024/1024*0.9}')
        ;;
    *)
        echo "Unsupported OS: $OS"
        exit 1
        ;;
esac
export MEM="${MEM}G"

gatk --java-options "-Xmx$MEM" SelectVariants -R "$REFERENCE" -V "$COHORT_VCF" -select-type SNP -O "$RAW_SNP" &
gatk --java-options "-Xmx$MEM" VariantFiltration -R "$REFERENCE" -V "$RAW_SNP" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "default_snp_filter" -O "$DEFAULT_SNP" &
gatk --java-options "-Xmx$MEM" VariantFiltration -R "$REFERENCE" -V "$DEFAULT_SNP" --filter-expression "AC < 4" --filter-name "frequency_filter" -O $FREQUENCY_SNP &
gatk --java-options "-Xmx$MEM" SelectVariants -R "$REFERENCE" -V "$FREQUENCY_SNP" --exclude-filtered -O "$FILTERED_SNP"

# rm $RAW_SNP $DEFAULT_SNP $FREQUENCY_SNP

# INDEL processing variables
RAW_INDEL="${PREFIX}.raw.INDELS.vcf.gz"
DEFAULT_INDEL="${PREFIX}.default.INDELS.vcf.gz"
FREQUENCY_INDEL="${PREFIX}.frequency.INDELS.vcf.gz"
FILTERED_INDEL="${PREFIX}.filtered.INDELS.vcf.gz"

mkfifo "$RAW_INDEL" "$DEFAULT_INDEL" "$FREQUENCY_INDEL"
# Process INDEL variants
gatk --java-options "-Xmx$MEM" SelectVariants -R "$REFERENCE" -V "$COHORT_VCF" -select-type INDEL -O "$RAW_INDEL" &
gatk --java-options "-Xmx$MEM" VariantFiltration -R "$REFERENCE" -V "$RAW_INDEL" --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "default_indel_filter" -O "$DEFAULT_INDEL" &
gatk --java-options "-Xmx$MEM" VariantFiltration -R "$REFERENCE" -V "$DEFAULT_INDEL" --filter-expression "AC < 4" --filter-name "frequency_filter" -O "$FREQUENCY_INDEL" &
gatk --java-options "-Xmx$MEM" SelectVariants -R "$REFERENCE" -V "$FREQUENCY_INDEL" --exclude-filtered -O "$FILTERED_INDEL"
# # Cleanup named pipes
# rm *.snps.tmp
