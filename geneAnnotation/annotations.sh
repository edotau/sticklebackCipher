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

# Set required paths, tools, and references
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
echo "
Creating chromSizes format...

    faSize -detailed $TARGET > $TCHROM
"
$faSize -detailed "$TARGET" > "$TCHROM"

# Step 2: Convert GTF to GenePred format
echo "
Converting $GTF to GenePred format...

    gtfToGenePred -genePredExt -includeVersion $GTF $GP_FILE
"
$gtfToGenePred -genePredExt -includeVersion "$GTF" "$GP_FILE"

# Step 3: Convert GenePred to PSL format
echo "
Converting $GP_FILE to PSL format...

    genePredToPsl $TCHROM $GP_FILE $PSL_FILE
"
$genePredToPsl "$TCHROM" "$GP_FILE" "$PSL_FILE"

# Step 4: Map PSL using chain file
echo "
Mapping $PSL_FILE using chain file...

    pslMap -chainMapFile -swapMap $PSL_FILE $CHAIN $MAPPED_PSL_FILE
"
$pslMap -chainMapFile -swapMap "$PSL_FILE" "$CHAIN" "$MAPPED_PSL_FILE"

# Step 5: Convert mapped PSL to BED format
echo "
Converting $MAPPED_PSL_FILE to BED format...

    pslToBed $MAPPED_PSL_FILE $OUTPUT_BED
"
$pslToBed "$MAPPED_PSL_FILE" "$OUTPUT_BED"

# Step 6: Convert BED to GenePred format
echo "
Converting $OUTPUT_BED to GenePred format...

    bedToGenePred $OUTPUT_BED $OUTPUT_GP
"
$bedToGenePred "$OUTPUT_BED" "$OUTPUT_GP"

# Step 7: Convert GenePred to GTF format
echo "
Converting $OUTPUT_GP to GTF format...

    genePredToGtf file $OUTPUT_GTF
"
$genePredToGtf file "$OUTPUT_GTF"

# Final Output
echo "Completed successfully!

Generated files:
    $OUTPUT_BED
    $OUTPUT_GP
    $OUTPUT_GTF
"

# Clean up intermediate files
rm -f "$TCHROM" "$GP_FILE" "$PSL_FILE"
