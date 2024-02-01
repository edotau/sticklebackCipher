#!/bin/bash
set -e
SCRIPT_NAME=$(basename "$0" | cut -d. -f1)

USAGE="
Usage:
    $SCRIPT_NAME <reference.fa> <read1.fq.gz> <read2.fq.gz>
"

REFERENCE=$1
READ_ONE=$2
READ_TWO=$3

if [ -z "$REFERENCE" ] || [ -z "$READ_ONE" ] || [ -z "$READ_TWO" ]; then
    echo "$USAGE"
    exit 1
fi

CURRENT_DIR=$(pwd)
PREFIX=$CURRENT_DIR/$(echo "$(basename "$(basename "$(basename \
    "$(basename "$READ_ONE" .fastq.gz)" \
    .fq.gz)" \
    .fastq)" \
    .fq)")

THREADS=$(nproc)
# Align reads to reference
TRIM_GALORE=/data/lowelab/edotau/software/TrimGalore-0.6.5/trim_galore

echo "Trimming reads with Trim Galore...
$TRIM_GALORE -o $CURRENT_DIR --basename $PREFIX --trim-n --max_n 0 --cores $(( $(nproc) > 4 ? 4 : $(nproc) )) --paired $READ_ONE $READ_TWO"

TRIMMED_FWD=${PREFIX}_val_1.fq.gz
TRIMMED_REV=${PREFIX}_val_2.fq.gz

FWD_READ=${PREFIX}_R1.fastq.gz
REV_READ=${PREFIX}_R2.fastq.gz

BBDUK="bbduk.sh"
# Second Trimming of adaptors:
echo "$BBDUK \
    in1=$TRIMMED_FWD in2=$TRIMMED_REV \
    out1=$FWD_READ out2=$REV_READ \
    minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
    ref=/data/lowelab/software/bbmap/resources/adapters.fa"

# MAPPING
echo -e "["$(date)"]\tAligning and sorting..."
SORT_THREADS=$(( $THREADS / 4 ))
BAM=${PREFIX}.bam

echo "
bwa mem -t $THREADS $REF $FWD_READ $REV_READ | samtools sort -@ $SORT_THREADS -T ${PREFIX}.tmp /dev/stdin > $BAM
"

echo "
samtools index $BAM
"