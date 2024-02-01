#!/bin/bash
set -eo pipefail

# Ensure tools are available
export PATH=$PATH:$GATK_PATH
export samtools=/usr/bin/samtools

# Validate arguments
if [ "$#" -ne 2 ]; then
    echo "
    Usage:
        ./$(basename "$0") <reference.fa> <sample.bam>
    "
    exit 1
fi

REFERENCE="$1"
BAM="$2"

PREFIX=$(basename "$BAM" .bam)

# Validate required tools
for cmd in gatk samtools; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in the current PATH."
        exit 1
    fi
done

PROCESSED_BAM="${PREFIX}.gatk.bam"
MARKED_DUPS=${PREFIX}.nodups.bam
THREADS=$(( $(nproc) > 12 ? 12 : $(nproc) ))

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

# Extract platform info
PL=$(samtools view "$BAM" | head -n 1 | cut -f 1 | cut -d ':' -f 3,4 | tr ':' '.')

# Pipeline processing
gatk --java-options "-Xmx$MEM" MarkDuplicates -I $BAM -O $MARKED_DUPS -M ${PREFIX}".dup.metrics" -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true

samtools view -b -@ "$THREADS" -f 1 $MARKED_DUPS | gatk --java-options "-Xmx$MEM" AddOrReplaceReadGroups -I /dev/stdin -O /dev/stdout -LB "WGS" -PL "$PL" -PU "ILUMINA" -SM "$PREFIX" | gatk --java-options "-Xmx$MEM" FixMateInformation -I /dev/stdin -O "$PROCESSED_BAM" -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true

# Index the final output BAM
samtools index "$PROCESSED_BAM"

# Variant calling
gatk --java-options "-Xmx$MEM" HaplotypeCaller -I "$PROCESSED_BAM" -O "${PREFIX}.g.vcf.gz" -R "$REFERENCE" -ERC GVCF --native-pair-hmm-threads "$THREADS"

# Assuming $gVCF should be ${PREFIX}.g.vcf.gz based on previous command
gatk --java-options "-Xmx$(echo $MEM | sed 's/G/g/')" GenotypeGVCFs -R "$REFERENCE" -V "${PREFIX}.g.vcf.gz" -O "${PREFIX}.vcf.gz"
