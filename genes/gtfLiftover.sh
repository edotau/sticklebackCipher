#!/bin/bash
set -e

# Input arguments
GTF=$1
TARGET_CHROM_SIZES=$2
CHAIN=$3
QUERY=$4

if [[ "$#" -lt 4 ]]; then
    echo "Usage:
        $0 annotations.gtf target.sizes target.chain query.fa
    "
    exit 1
fi

# Set required paths, tools and references
KENT_UTILS="/home"
gtfToGenePred="$KENT_UTILS/gtfToGenePred"
genePredToPsl="$KENT_UTILS/genePredToPsl"
pslMap="$KENT_UTILS/pslMap"
pslToBed="$KENT_UTILS/pslToBed"
bedToGenePred="$KENT_UTILS/bedToGenePred"
genePredToGtf="$KENT_UTILS/genePredToGtf"

# Check if all tools exist and are executable
if ! [ -x "$gtfToGenePred" ]; then
    echo "Error: gtfToGenePred not found or not executable at $gtfToGenePred"
    exit 1
fi

if ! [ -x "$genePredToPsl" ]; then
    echo "Error: genePredToPsl not found or not executable at $genePredToPsl"
    exit 1
fi

if ! [ -x "$pslMap" ]; then
    echo "Error: pslMap not found or not executable at $pslMap"
    exit 1
fi

if ! [ -x "$pslToBed" ]; then
    echo "Error: pslToBed not found or not executable at $pslToBed"
    exit 1
fi

if ! [ -x "$bedToGenePred" ]; then
    echo "Error: bedToGenePred not found or not executable at $bedToGenePred"
    exit 1
fi

if ! [ -x "$genePredToGtf" ]; then
    echo "Error: genePredToGtf not found or not executable at $genePredToGtf"
    exit 1
fi

# Paths to reference files

# Check if the chrom sizes file exists
if ! [ -f "$TARGET_CHROM_SIZES" ]; then
    echo "Error: Chrom sizes file not found at $TARGET_CHROM_SIZES"
    exit 1
fi

# Derived variables
PREFIX=$(basename "$GTF" .gtf)
OUTPUT=$(basename "$QUERY" .fa)
GP_FILE="${PREFIX}.gp"
PSL_FILE="${PREFIX}.psl"
MAPPED_PSL_FILE="${OUTPUT}.psl"
OUTPUT_BED="${OUTPUT}.bed"
OUTPUT_GP="${OUTPUT}.gp"
OUTPUT_GTF="${OUTPUT}.gtf"

# Step 1: Convert GTF to GenePred format
if ! [ -e "$GP_FILE" ]; then
    echo "Converting $GTF to GenePred format..."
    $gtfToGenePred -genePredExt -includeVersion "$GTF" "$GP_FILE"
fi

# Step 2: Convert GenePred to PSL format
if ! [ -e "$PSL_FILE" ]; then
    echo "Converting $GP_FILE to PSL format..."
    $genePredToPsl "$TARGET_CHROM_SIZES" "$GP_FILE" "$PSL_FILE"
fi

# Step 3: Map PSL using chain file
if ! [ -e "$MAPPED_PSL_FILE" ]; then
    echo "Mapping $PSL_FILE using chain file..."
    $pslMap -chainMapFile -swapMap "$PSL_FILE" "$CHAIN" "$MAPPED_PSL_FILE"
fi

# Step 4: Convert mapped PSL to BED format
if ! [ -e "$OUTPUT_BED" ]; then
    echo "Converting $MAPPED_PSL_FILE to BED format..."
    $pslToBed "$MAPPED_PSL_FILE" "$OUTPUT_BED"
fi

# Step 5: Convert BED to GenePred format
if ! [ -e "$OUTPUT_GP" ]; then
    echo "Converting $OUTPUT_BED to GenePred format..."
    $bedToGenePred "$OUTPUT_BED" "$OUTPUT_GP"
fi

# Step 6: Convert GenePred to GTF format
if ! [ -e "$OUTPUT_GTF" ]; then
    echo "Converting $OUTPUT_GP to GTF format..."
    $genePredToGtf file "$OUTPUT_GTF"
fi

# Final Output
echo "All operations completed successfully."
echo "Generated files:"
echo "$OUTPUT_BED"
echo "$OUTPUT_GP"
echo "$OUTPUT_GTF"

# Clean up intermediate files
rm -f "$GP_FILE" "$PSL_FILE"
