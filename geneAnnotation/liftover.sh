#!/bin/bash -e

# Input arguments
TARGET=$1
QUERY=$2
TCHAIN=$3

if [[ "$#" -lt 3 ]]; then
    echo "Usage:
        $(basename "$0") target.fa query.fa target.chain
    "
    exit 1
fi

# Set required paths, tools, and references
KENT_UTILS="/home/kentUtils/bin"
HTSLIB_UTILS="/home/htslib/bin"

# Append $KENT_UTILS and $HTSLIB_UTILS or set binary directly
faToTwoBit="$KENT_UTILS/faToTwoBit"
faSize="$KENT_UTILS/faSize"
chainStitchId="$KENT_UTILS/chainStitchId"
chainSwap="$KENT_UTILS/chainSwap"
chainSort="$KENT_UTILS/chainSort"
chainPreNet="$KENT_UTILS/chainPreNet"
chainNet="$KENT_UTILS/chainNet"
netSyntenic="$KENT_UTILS/netSyntenic"
netChainSubset="$KENT_UTILS/netChainSubset"
chainToPsl="$KENT_UTILS/chainToPsl"
netToAxt="$KENT_UTILS/netToAxt"
axtSort="$KENT_UTILS/axtSort"
samtools="$HTSLIB_UTILS/samtools"

# Check if all tools exist and are executable
for i in "$faToTwoBit" "$faSize" "$chainStitchId" "$chainSwap" "$chainSort" "$chainPreNet" "$chainNet" "$netSyntenic" "$netChainSubset" "$chainToPsl" "$netToAxt" "$axtSort" "$samtools"; do
    if ! [ -x "$i" ]; then
        echo "Error: $(basename "$i") not found or not executable at $i"
        exit 1
    fi
done

# Derived variables
TNAME=$(basename "$TARGET" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')
QNAME=$(basename "$QUERY" | sed -E 's/\.(fa|fa\.gz|fasta|fasta\.gz)$//')

BITTARGET="${TNAME}.2bit"
BITQUERY="${QNAME}.2bit"
TCHROM="${TNAME}.sizes"
QCHROM="${QNAME}.sizes"

# Step 1: Convert target and query to 2bit format if not already done
echo "
Converting $TARGET and $QUERY to 2bit format ...

    $faToTwoBit "$TARGET" "$BITTARGET"
    $faToTwoBit "$QUERY" "$BITQUERY"
"
$faToTwoBit "$TARGET" "$BITTARGET"
$faToTwoBit "$QUERY" "$BITQUERY"

# Step 2: Create chromSizes for target and query if not already done
echo "
Creating chromSizes format...

    $faSize -detailed "$TARGET" > "$TCHROM"
    $faSize -detailed "$QUERY" > "$QCHROM"
"
$faSize -detailed "$TARGET" > "$TCHROM"
$faSize -detailed "$QUERY" > "$QCHROM"

# Step 3: Swap target best chains to be query-referenced
QRY_TGT="${QNAME}.${TNAME}.tBest.chain"
echo "
Swapping chains to query-referenced ...

    chainStitchId $TCHAIN /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin ${QNAME}_${TNAME}.tBest.chain
"
$chainStitchId "$TCHAIN" /dev/stdout | $chainSwap /dev/stdin /dev/stdout | $chainSort /dev/stdin "$QRY_TGT"

# Step 4: Net query-referenced chains to get reciprocal best net
QRY_TGT_NET="${QNAME}.${TNAME}.rBest.net.gz"
echo "
Generating reciprocal best net for query-referenced chains...

    chainPreNet $QRY_TGT $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > $QRY_TGT_NET
"
$chainPreNet "$QRY_TGT" "$QCHROM" "$TCHROM" /dev/stdout | \
$chainNet -minSpace=1 -minScore=0 /dev/stdin "$QCHROM" "$TCHROM" /dev/stdout /dev/null | \
$netSyntenic /dev/stdin /dev/stdout | gzip -c > "$QRY_TGT_NET"

# Step 5: Extract query-referenced reciprocal best chain
QRY_TGT_NETCHAIN="${QNAME}.${TNAME}.rBest.netted.chain.gz"
echo "
Extracting reciprocal best chain for query-referenced chains...

    netChainSubset $QRY_TGT_NET $QRY_TGT /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > $QRY_TGT_NETCHAIN
"
$netChainSubset "$QRY_TGT_NET" "$QRY_TGT" /dev/stdout | $chainStitchId /dev/stdin /dev/stdout | gzip -c > "$QRY_TGT_NETCHAIN"

# Step 6: Swap to get target-referenced reciprocal best chain
TGT_QRY_SWAP_CHAIN="${TNAME}.${QNAME}.rBest.netted.chain.gz"
echo "
Swapping to target-referenced reciprocal best chain...

    chainSwap $QRY_TGT_NETCHAIN /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_CHAIN
"
$chainSwap "$QRY_TGT_NETCHAIN" /dev/stdout | $chainSort /dev/stdin /dev/stdout | gzip -c > "$TGT_QRY_SWAP_CHAIN"

# Step 7: Net target-referenced chains to get target-referenced reciprocal best net
TGT_QRY_SWAP_NET="${TNAME}.${QNAME}.rbest.net.gz"
echo "
Generating target-referenced reciprocal best net...

    chainPreNet $TGT_QRY_SWAP_CHAIN $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_NET
"
$chainPreNet "$TGT_QRY_SWAP_CHAIN" "$TCHROM" "$QCHROM" /dev/stdout | \
$chainNet -minSpace=1 -minScore=0 /dev/stdin "$TCHROM" "$QCHROM" /dev/stdout /dev/null | \
$netSyntenic /dev/stdin /dev/stdout | gzip -c > "$TGT_QRY_SWAP_NET"

# Clean up intermediate files
rm -f "$QRY_TGT"

# Step 8: Generate AXT and SAM/BAM files from nets
TGT_QRY_SWAP_NET_AXT="${TNAME}.${QNAME}.rBest.netted.chain.axt.gz"
echo "
Converting nets to AXT format...

    netToAxt $TGT_QRY_SWAP_NET $TGT_QRY_SWAP_CHAIN $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_NET_AXT
"
$netToAxt "$TGT_QRY_SWAP_NET" "$TGT_QRY_SWAP_CHAIN" "$BITTARGET" "$BITQUERY" /dev/stdout | $axtSort /dev/stdin /dev/stdout | gzip -c > "$TGT_QRY_SWAP_NET_AXT"

rm -f "$BITTARGET" "$BITQUERY" "$TCHROM" "$QCHROM" "$QRY_TGT" "$QRY_TGT_NET" "$QRY_TGT_NETCHAIN" "$TGT_QRY_SWAP_CHAIN"

# Success message
echo "Completed successfully!

Generated files:
    $QRY_TGT_NET
    $QRY_TGT_NETCHAIN
    $TGT_QRY_SWAP_CHAIN
    $TGT_QRY_SWAP_NET
    $TGT_QRY_SWAP_NET_AXT
"