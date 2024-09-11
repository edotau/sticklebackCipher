#!/bin/bash -e

target=$1
query=$2
DIR=${3%/} # Folder containing chains to merge

if [[ "$#" -lt 3 ]]; then
    echo "Usage:
        $(basename "$0") target.fa query.fa dir/chain/
    "
    exit 1
fi

# Set required paths, tools, and references
KENT_UTILS="/Users/edotau/src/github.com/kentUtils/bin"

# Append $KENT_UTILS or set binary directly
faToTwoBit="$KENT_UTILS/faToTwoBit"
twoBitInfo="$KENT_UTILS/twoBitInfo"
chainMergeSort="$KENT_UTILS/chainMergeSort"
chainSort="$KENT_UTILS/chainSort"
chainSplit="$KENT_UTILS/chainSplit"

TNAME=$(basename "$target" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')
TBIT=${TNAME}.2bit
TCROM=${TNAME}.sizes

QNAME=$(basename "$query" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')
QCHROM=${QNAME}.sizes
QBIT=${QNAME}.2bit

OUTPUT="${TNAME}.${QNAME}.all.sorted.chain"

# Generate 2bit files if they dont exist
if [[ ! -f "$TBIT" ]]; then
    echo "Generating target.2bit...
    
        faToTwoBit ""$target"" ""$TBIT""
    "
    $faToTwoBit "$target" "$TBIT"
fi
if [[ ! -f "$QBIT" ]]; then
    echo "Generating query chrom sizes...

        faToTwoBit $query $QBIT
    "
    faToTwoBit "$query" "$QBIT"
fi

# Generate chrom sizes if they don't exist
if [[ ! -f "$TCROM" ]]; then
    echo "Generating target chrom sizes...
    
        $twoBitInfo $TBIT $TCHROM
    "
    $twoBitInfo "$TBIT" "$TCROM"
fi
if [[ ! -f "$QCHROM" ]]; then
    echo "Generating query chrom sizes...

        $twoBitInfo $QBIT $QCHROM
    "
    $twoBitInfo "$QBIT" "$QCHROM"
fi

# Merge and sort chains
echo "
Merging and splitting chains...

    chainMergeSort $DIR/*.chain | $chainSplit chainMerge stdin -lump=50
"
$chainMergeSort "$DIR"/*.chain | $chainSplit chainMerge stdin -lump=50

echo "
Sorting merged chains...

    cat chainMerge/*.chain | $chainSort /dev/stdin $OUTPUT
"
cat chainMerge/*.chain | $chainSort /dev/stdin "$OUTPUT"

#Clean up intermediate files
rm "$TCROM" "$QCHROM" "$TBIT" "$QBIT"
rm -rf chainMerge

echo "Completed successfully!

Output written to $OUTPUT !
"
