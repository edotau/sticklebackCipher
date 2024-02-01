#!/bin/bash
set -e pipefail

REFERENCE=$1
BAM=$2

USAGE="
Usage:
    $(basename "$0" | cut -d. -f1) <reference.fa> <sample.bam>
"

if [ -z "$REFERENCE" ] || [ -z "$BAM" ]; then
    echo "$USAGE"
    exit 1
fi

PREFIX=$(pwd)/$(basename "$BAM" .bam)

# Define required arguments
LIB="WGS"
PU="ILUMINA"
ID=$PREFIX

# Extract PL information from BAM header
PL=$(samtools view "$BAM" | head -n 1 | cut -f 1 | cut -d ':' -f 3,4 | tr ':' '.')



# Pipeline for MarkDuplicates, AddOrReplaceReadGroups, and FixMateInformation
echo "
gatk --java-options "-Xmx16G" MarkDuplicates -I "$BAM" -O /dev/stdout -M /dev/null -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true |
samtools view -b -@ $(( $(nproc) > 4 ? 4 : $(nproc) )) -f 1 - |
gatk --java-options "-Xmx16G" AddOrReplaceReadGroups -I /dev/stdin -O /dev/stdout -LB "$LIB" -PL "$PL" -PU "$PU" -SM "$ID" |
gatk --java-options "-Xmx16G" FixMateInformation -I /dev/stdin -O "${PREFIX}.gatk.bam" -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR "$DIR/fixmate.tmp"
"

# Index the final output BAM
echo "
samtools index "$DIR/${PREFIX}.gatk.valid.bam"
"

echo "Ready to run GATK on processed $PREFIX BAM"
