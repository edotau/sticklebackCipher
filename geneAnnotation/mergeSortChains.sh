#!/bin/bash -e

targetTwoBit=$1
queryTwoBit=$2
DIR=${3%/} # Folder containing chains to merge

if [[ "$#" -lt 3 ]]; then
    echo "Usage:
        $(basename "$0") target.2bit query.2bit chains/
    "
    exit 1
fi

# Set required paths, tools, and references
KENT_UTILS="/home/kentUtils/bin"

# Append $KENT_UTILS or set binary directly
twoBitInfo="$KENT_UTILS/twoBitInfo"
chainMergeSort="$KENT_UTILS/chainMergeSort"
chainSort="$KENT_UTILS/chainSort"
chainSplit="$KENT_UTILS/chainSplit"

TNAME=$(basename ${targetTwoBit} .2bit)
TCROM=${TNAME}.sizes
QNAME=$(basename ${queryTwoBit} .2bit)
QCHROM=${QNAME}.sizes
OUTPUT="${TNAME}.${QNAME}.all.sorted.chain"

# Generate chrom sizes if they don't exist
if [[ ! -f "$TCROM" ]]; then
    echo "Generating target chrom sizes...
    
        $twoBitInfo $targetTwoBit $TCHROM
    "
    $twoBitInfo "$targetTwoBit" "$TCROM"
fi

if [[ ! -f "$QCHROM" ]]; then
    echo "Generating query chrom sizes...

        $twoBitInfo $queryTwoBit $QCHROM
    "
    $twoBitInfo "$queryTwoBit" "$QCHROM"
fi

# Merge and sort chains
echo "
Merging and splitting chains...

    $chainMergeSort $DIR/*.chain | $chainSplit chainMerge stdin -lump=50
"
$chainMergeSort $DIR/*.chain | $chainSplit chainMerge stdin -lump=50

echo "
Sorting merged chains...

    cat chainMerge/*.chain | $chainSort /dev/stdin $OUTPUT
"
cat chainMerge/*.chain | $chainSort /dev/stdin $OUTPUT

echo "Completed successfully!
Output written to $OUTPUT.
"
