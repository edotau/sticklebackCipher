#!/bin/sh
module load gcc
hgLoadChain=/data/lowelab/cl454/bin/x86_64/hgLoadChain
bedToBigBed=/data/lowelab/cl454/bin/x86_64/bedToBigBed
chain=$1
chromSizes=$2
PREFIX=$(basename $chain .chain)

$hgLoadChain -noBin -test hg38 bigChain $chain
sed 's/.000000//' chain.tab | awk 'BEGIN {OFS="\t"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}' > ${PREFIX}.bigChain
$bedToBigBed -type=bed6+6 -as=/data/lowelab/edotau/scripts/bigChain.as -tab ${PREFIX}.bigChain $chromSizes ${PREFIX}.bb
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5, $4}' link.tab | sort -k1,1 -k2,2n > ${PREFIX}.bigChain.bigLink
$bedToBigBed -type=bed4+1 -as=/data/lowelab/edotau/scripts/bigLink.as -tab ${PREFIX}.bigChain.bigLink $chromSizes ${PREFIX}.link.bb
