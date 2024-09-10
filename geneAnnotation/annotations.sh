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
bedToGenePred="$KENT_UTILS/bedToGenePred"
bedToPsl="$KENT_UTILS/bedToPsl"
gtfToGenePred="$KENT_UTILS/gtfToGenePred"
genePredToGtf="$KENT_UTILS/genePredToGtf"
genePredToBed="$KENT_UTILS/genePredToBed"
faSize="$KENT_UTILS/faSize"
pslMap="$KENT_UTILS/pslMap"
pslToBed="$KENT_UTILS/pslToBed"

# Check if all tools exist and are executable
for i in "$gtfToGenePred" "$bedToPsl" "$pslMap" "$pslToBed" "$bedToGenePred" "$genePredToGtf" "$genePredToBed" "$faSize"; do
    if ! [ -x "$i" ]; then
        echo "Error: $(basename "$i") not found or not executable at $i"
        exit 1
    fi
done

# Derived variables
PREFIX=$(basename "$GTF" .gtf)
TCHROM=$(basename "$TARGET" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//').sizes
GP_FILE="${PREFIX}.gp"
PSL_FILE="${PREFIX}.psl"
MAPPED_PSL_FILE="${OUTPUT}.psl"
OUTPUT_BED="${OUTPUT}.bed"
OUTPUT_GP="${OUTPUT}.gp"
OUTPUT_GTF="${OUTPUT}.gtf"

# Step 1: Create target chromSizes from target .fa
if [[ ! -f "$TCROM" ]]; then
    echo "Generating target chrom sizes...
    
        $faSize -detailed $TARGET > $TCHROM
    "
    $faSize -detailed "$TARGET" > "$TCHROM"
fi

# Step 2: Convert GTF to GenePred format
echo "
Converting $GTF to PSL format...

    $gtfToGenePred -genePredExt -includeVersion $GTF /dev/stdout | $genePredToBed /dev/stdin  /dev/stdout | $bedToPsl -tabs -keepQuery $TCHROM /dev/stdin $PSL_FILE
"
$gtfToGenePred -genePredExt -includeVersion "$GTF" /dev/stdout | $genePredToBed /dev/stdin /dev/stdout | $bedToPsl -tabs -keepQuery "$TCHROM" /dev/stdin "$PSL_FILE"


# Step 3: Map PSL using chain file
echo "
Mapping $PSL_FILE using chain file...

    $pslMap -chainMapFile -swapMap $PSL_FILE $CHAIN $MAPPED_PSL_FILE
"
$pslMap -chainMapFile -swapMap "$PSL_FILE" "$CHAIN" "$MAPPED_PSL_FILE"

# Step 4: Convert mapped PSL to BED format
echo "
Converting $MAPPED_PSL_FILE to GTF format...

    $pslToBed $MAPPED_PSL_FILE /dev/stdout | $bedToGenePred /dev/stdin /dev/stdout | $genePredToGtf file /dev/stdin $OUTPUT_GTF
"
$pslToBed "$MAPPED_PSL_FILE" /dev/stdout | $bedToGenePred /dev/stdin /dev/stdout | $genePredToGtf file /dev/stdin "$OUTPUT_GTF"

# Final Output
echo "Completed successfully!

Generated files:
    $OUTPUT_BED
    $OUTPUT_GP
    $OUTPUT_GTF
"

