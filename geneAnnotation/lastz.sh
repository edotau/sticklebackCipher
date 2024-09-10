#!/bin/bash -e

TARGET=$1 # Target/reference fa MUST be a single record fasta sequence
QUERY=$2  # Query fa can be multi records

if [[ "$#" -lt 2 ]]; then
    echo "Usage:
        $(basename "$0") target.fa query.fa
    "
    exit 0
fi

# Set lastz or add to path
KENT_UTILS="/home/kentUtils/bin"

lastz="/usr/local/bin/lastz"
axtChain="$KENT_UTILS/axtChain"
chainSort="$KENT_UTILS/chainSort"

# Check if all tools exist and are executable
for i in "$lastz" "$axtChain"; do
    if ! [ -x "$i" ]; then
        echo "Error: $(basename "$i") not found or not executable at $i"
        exit 1
    fi
done

# Determine Scoring matrix
MATRIX="human-chimp-matrix.txt"
touch $MATRIX && rm -f $MATRIX
echo "A     C     G     T
    A   91  -114   -31  -123
    C -114   100  -125   -31
    G  -31  -125   100  -114
    T -123   -31  -114    91
" > $MATRIX

PREFIX=$(basename "$TARGET" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')_$(basename "$QUERY" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')

# Output files
axt=${PREFIX}.axt
chainFile=${PREFIX}.chain

# Step 1: Run lastz alignment
echo "
Running lastz alignment...

    $lastz $TARGET $QUERY --format=axt --output=$axt --scores=$MATRIX O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
"
$lastz "$TARGET" "$QUERY" --format=axt --output="$axt" --scores="$MATRIX" O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
echo "Finished lastz alignment!"

# Step 2: Convert AXT to chain format
echo "
Converting AXT to chain format...

    $axtChain -linearGap=medium -scoreScheme=$MATRIX $axt -faT $TARGET -faQ $QUERY /dev/stdout | $chainSort /dev/stdin $chainFile
"
$axtChain -linearGap=medium -scoreScheme="$MATRIX" "$axt" -faT "$TARGET" -faQ "$QUERY" /dev/stdout | $chainSort /dev/stdin "$chainFile"
echo "Finished converting axt to chain!"

# Clean up intermediate files
echo "Cleaning up files..."
rm -f $MATRIX

# Success message
echo "Completed successfully!

Generated files:
    $axt
    $chainFile
"
