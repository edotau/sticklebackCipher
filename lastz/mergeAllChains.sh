#!/bin/sh
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --time=0-01:00:00
#Requires target.2bit and query.2bit and and directory to the chain files
module add kentUtils/v302-gcb01
targetTwoBit=$1
queryTwoBit=$2
DIR=${3%/}
PREFIX=$(basename $targetTwoBit .2bit)to$(basename ${queryTwoBit^} .2bit)
target=$(basename $targetTwoBit .2bit).sizes
query=$(basename $queryTwoBit .2bit).sizes

twoBitInfo $targetTwoBit $target
twoBitInfo $queryTwoBit $query

chainMergeSort $DIR/*.chain | chainSplit chainMerge stdin -lump=50
cat chainMerge/*.chain > ${PREFIX}.all.chain
chainSort ${PREFIX}.all.chain ${PREFIX}.all.sorted.chain
rm -r chainMerge

chainPreNet ${PREFIX}.all.sorted.chain $target $query ${PREFIX}ChainPreNet.chain
net=${PREFIX}.net
chainNet ${PREFIX}ChainPreNet.chain $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin $net
${PREFIX}_syntenic_netted_chain.net
netFilter -minScore=500000 $net > ${PREFIX}Filtered.net

netChainSubset ${PREFIX}Filtered.net ${PREFIX}.all.sorted.chain ${PREFIX}Netted.chain
axt=${PREFIX}NettedChain.axt
netToAxt ${PREFIX}Filtered.net ${PREFIX}.all.sorted.chain $bitTarget $bitQuery $axt
rm $target $query ${PREFIX}ChainPreNet.chain
