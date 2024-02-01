#!/bin/bash
set -eo pipefail

# Ensure tools are available
export PATH=$PATH:$GATK_PATH
export samtools=/usr/bin/samtools

if [ "$#" -ne 2 ]; then
    echo "Usage: $(basename "$0") <reference.fa> <directory>"
    exit 1
fi

REFERENCE="$1"
DIR="${2%/}"

# Extract PREFIX by removing potential reference file extensions
PREFIX=$(basename "$REFERENCE" | sed -E 's/\.(fa|fasta)(\.gz)?$//')

# Verify required commands are in PATH
for cmd in gatk samtools; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in the current PATH."
        exit 1
    fi
done

SAMPLES="${DIR}/${PREFIX}.list"
: > "$SAMPLES"

for i in "$DIR"/*.g.vcf.gz; do
    echo "$i" >> "$SAMPLES"
done

MERGED="${DIR}/${PREFIX}.merged.g.vcf.gz"
VCF="${DIR}/${PREFIX}.cohort.genotyped.vcf.gz"

# Combine GVCFs
gatk --java-options "-Xmx32g" CombineGVCFs -R "$REFERENCE" --variant "$SAMPLES" -O "$MERGED"

# Genotype GVCFs
gatk --java-options "-Xmx32g" GenotypeGVCFs -R "$REFERENCE" -V "$MERGED" -O "$VCF"
