#!/bin/bash -e

# Input arguments
TARGET=$1
CHAIN=$2
GTF=$3
OUTPUT=$4

if [[ "$#" -lt 4 ]]; then
    echo "Usage:
        $(basename "$0") target.fa target.chain annotations.gtf output
    "
    exit 1
fi

# Set required paths, tools and references
KENT_UTILS="/home/kentUtils/bin"

# Append $KENT_UTILS or set binary directly
gtfToGenePred="$KENT_UTILS/gtfToGenePred"
genePredToPsl="$KENT_UTILS/genePredToPsl"
pslMap="$KENT_UTILS/pslMap"
pslToBed="$KENT_UTILS/pslToBed"
bedToGenePred="$KENT_UTILS/bedToGenePred"
genePredToGtf="$KENT_UTILS/genePredToGtf"
faSize="$KENT_UTILS/faSize"

# Check if all tools exist and are executable
for i in "$gtfToGenePred" "$genePredToPsl" "$pslMap" "$pslToBed" "$bedToGenePred" "$genePredToGtf" "$faSize"; do
    if ! [ -x "$i" ]; then
        echo "Error: $(basename "$i") not found or not executable at $i"
        exit 1
    fi
done

# Derived variables
PREFIX=$(basename "$GTF" .gtf)
TCHROM=$(basename "$TARGET" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')
GP_FILE="${PREFIX}.gp"
PSL_FILE="${PREFIX}.psl"
MAPPED_PSL_FILE="${OUTPUT}.psl"
OUTPUT_BED="${OUTPUT}.bed"
OUTPUT_GP="${OUTPUT}.gp"
OUTPUT_GTF="${OUTPUT}.gtf"

# Step 1: Create target chromSizes from target .fa
if ! [ -e "$TCHROM" ]; then
    echo "Creating chromSizes format..."
    $faSize -detailed "$TARGET" > "$TCHROM"
fi

# Step 2: Convert GTF to GenePred format
if ! [ -e "$GP_FILE" ]; then
    echo "Converting $GTF to GenePred format..."
    $gtfToGenePred -genePredExt -includeVersion "$GTF" "$GP_FILE"
fi

# Step 3: Convert GenePred to PSL format
if ! [ -e "$PSL_FILE" ]; then
    echo "Converting $GP_FILE to PSL format..."
    $genePredToPsl "$TCHROM" "$GP_FILE" "$PSL_FILE"
fi

# Step 4: Map PSL using chain file
if ! [ -e "$MAPPED_PSL_FILE" ]; then
    echo "Mapping $PSL_FILE using chain file..."
    $pslMap -chainMapFile -swapMap "$PSL_FILE" "$CHAIN" "$MAPPED_PSL_FILE"
fi

# Step 5: Convert mapped PSL to BED format
if ! [ -e "$OUTPUT_BED" ]; then
    echo "Converting $MAPPED_PSL_FILE to BED format..."
    $pslToBed "$MAPPED_PSL_FILE" "$OUTPUT_BED"
fi

# Step 6: Convert BED to GenePred format
if ! [ -e "$OUTPUT_GP" ]; then
    echo "Converting $OUTPUT_BED to GenePred format..."
    $bedToGenePred "$OUTPUT_BED" "$OUTPUT_GP"
fi

# Step 7: Convert GenePred to GTF format
if ! [ -e "$OUTPUT_GTF" ]; then
    echo "Converting $OUTPUT_GP to GTF format..."
    $genePredToGtf file "$OUTPUT_GTF"
fi

# Final Output
echo "Completed successfully!

Generated files:
    $OUTPUT_BED
    $OUTPUT_GP
    $OUTPUT_GTF
"

# Clean up intermediate files
rm -f "$TCHROM" "$GP_FILE" "$PSL_FILE"
