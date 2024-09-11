#!/bin/bash -e

TARGET=$1 # Target/Reference fasta MUST be a single sequence record
QUERY=$2  # Query fa can be single or multi records

if [[ "$#" -lt 2 ]]; then
    echo "
    Usage:
        $(basename "$0") target.fa query.fa
    "
    exit 0
fi

# Set lastz or add to path
KENT_UTILS="/Users/edotau/src/github.com/kentUtils/bin"

lastz="/usr/local/bin/lastz"
axtChain="$KENT_UTILS/axtChain"
chainSort="$KENT_UTILS/chainSort"

# Check if all tools exist and are executable
for i in "$lastz" "$axtChain" "$chainSort"; do
    if ! [ -x "$i" ]; then
        echo "Error: $(basename "$i") not found or not executable at $i"
        exit 1
    fi
done

PREFIX=$(basename "$TARGET" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')_$(basename "$QUERY" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')

# Determine Scoring matrix
MATRIX="scores.mat"
if ! [ -e "$MATRIX" ]; then
    echo "\
hsp_threshold      = 3000
gapped_threshold   = 3000
x_drop             = 900
y_drop             = 45600
gap_open_penalty   = 600
gap_extend_penalty = 150

    A    C    G    T
A   90 -330 -236 -356
C -330  100 -318 -236
G -236 -318  100 -330
T -356 -236 -330   90" > $MATRIX
fi

# Output files
AXT=${PREFIX}.axt
chainFile=${PREFIX}.chain

# Step 1: Run lastz alignment and convert axt format to chain
echo "
Starting lastz alignment...

    lastz $TARGET $QUERY --format=axt --output=$AXT--scores=$MATRIX O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
"
$lastz "$TARGET" "$QUERY" --format=axt  --scores="$MATRIX" --output="$AXT" O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
echo "Finished lastz alignment!"

# Step 2: Convert AXT to chain format
echo "
Converting AXT to chain format and perform a sort for further processing...

    axtChain -linearGap=medium -scoreScheme=$MATRIX $AXT -faT $TARGET -faQ $QUERY /dev/stdout | chainSort /dev/stdin $chainFile
"
$axtChain -linearGap=medium -scoreScheme="$MATRIX" "$AXT" -faT "$TARGET" -faQ "$QUERY" /dev/stdout | $chainSort /dev/stdin "$chainFile"

# Clean up intermediate files
rm -f $MATRIX "$AXT"

# Success message
echo "
Completed successfully!

Generated file:
    $chainFile
"
