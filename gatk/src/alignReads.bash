#!/bin/bash
set -eo pipefail

# Ensure tools are available
export PATH=$PATH:$BBDUK_PATH:$TRIM_GALORE_PATH
export ADAPTERS=$BBDUK_PATH/resources/adapters.fa

export samtools=/usr/bin/samtools
export bwa=/usr/bin/bwa
export bbduk="$BBDUK_PATH/bbduk.sh"

# Validate the number of arguments
if [ "$#" -ne 3 ]; then
    echo "
    Usage:
        ./$(basename "$0") <reference.fa> <read1.fastq.gz> <read2.fastq.gz>
    "
    exit 1
fi

REFERENCE="$1"
READ_ONE="$2"
READ_TWO="$3"

# Extract the base prefix, handling different extensions properly
PREFIX=$(basename "$READ_ONE" | sed -E 's/_R[12]//; s/\.(fastq|fq)(\.gz)?$//')

# Define the output directory based on READ_ONE path
OUTPUT_DIR=$(dirname "$READ_ONE")

for cmd in trim_galore $bbduk samtools; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in the PATH."
        exit 1
    fi
done

# Set the number of processing threads, with a maximum of 4
THREADS=$(nproc)

# Perform read trimming with Trim Galore
echo "["$(date)"]\tTrimming reads with Trim Galore..."
trim_galore --output_dir "$OUTPUT_DIR" --basename "$PREFIX" --trim-n --max_n 0 --cores $(( THREADS > 4 ? 4 : THREADS )) --paired "$READ_ONE" "$READ_TWO"

# Define paths for Trim Galore output
TRIMMED_FWD="${OUTPUT_DIR}/${PREFIX}_val_1.fq.gz"
TRIMMED_REV="${OUTPUT_DIR}/${PREFIX}_val_2.fq.gz"

# Adapter trimming and quality filtering with BBDuk
echo "["$(date)"]\tRemoving adapters with BBDuk..."

$bbduk \
    in1="$TRIMMED_FWD" in2="$TRIMMED_REV" \
    out1="${OUTPUT_DIR}/${PREFIX}_R1.fastq.gz" out2="${OUTPUT_DIR}/${PREFIX}_R2.fastq.gz" \
    minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
    ref="$ADAPTERS"

# MAPPING with bwa and sorting with samtools
echo -e "["$(date)"]\tAligning and sorting..."

BAM="${OUTPUT_DIR}/${PREFIX}.bam"
SORT_THREADS=$(( $THREADS / 2 )) # Adjust based on your system's capabilities

bwa mem -t "$THREADS" "$REFERENCE" "${OUTPUT_DIR}/${PREFIX}_R1.fastq.gz" "${OUTPUT_DIR}/${PREFIX}_R2.fastq.gz" | \
samtools sort -@ "$SORT_THREADS" -T "${OUTPUT_DIR}/${PREFIX}.tmp" -o "$BAM"
samtools index "$BAM"
