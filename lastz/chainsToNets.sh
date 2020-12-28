#!/bin/sh
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --time=0-02
#SBATCH --nodes=1
set -e

###INPUTS:
targetContigs=$1
queryContigs=$2
DIR=${3%/}
PREFIX=$4

export PATH=/data/lowelab/edotau/kentUtils/:/data/lowelab/edotau/bin/envs/htslib/bin/:$PATH

targetPREFIX=$(basename $targetContigs .fa)
queryPREFIX=$(basename $queryContigs .fa)

target=${targetPREFIX}.sizes
query=${queryPREFIX}.sizes

bitTarget=${targetPREFIX}.2bit
bitQuery=${queryPREFIX}.2bit

net=${PREFIX}.net
axt=${PREFIX}.netted.chains.axt

echo "
faSize -detailed $targetContigs > $target"
faSize -detailed $targetContigs > $target

echo "
faSize -detailed $queryContigs > $query"
faSize -detailed $queryContigs > $query

echo "
faToTwoBit $targetContigs $bitTarget"
faToTwoBit $targetContigs $bitTarget

echo "
faToTwoBit $queryContigs $bitQuery"
faToTwoBit $queryContigs $bitQuery

echo "
chainMergeSort $DIR/*.chain | chainSplit chainMerge stdin -lump=50"
chainMergeSort $DIR/*.chain | chainSplit chainMerge stdin -lump=50

echo "
cat chainMerge/*.chain > ${PREFIX}.all.chain"
cat chainMerge/*.chain > ${PREFIX}.all.chain

echo "
chainSort ${PREFIX}.all.chain ${PREFIX}.all.sorted.chain; cat ${PREFIX}.all.sorted.chain | grep -v '#' > ${PREFIX}.all.chain; rm ${PREFIX}.all.sorted.chain"
chainSort ${PREFIX}.all.chain ${PREFIX}.all.sorted.chain; cat ${PREFIX}.all.sorted.chain | grep -v '#' > ${PREFIX}.all.chain; rm ${PREFIX}.all.sorted.chain
rm -r chainMerge

echo "
chainPreNet ${PREFIX}.all.chain $target $query ${PREFIX}ChainPreNet.chain"
chainPreNet ${PREFIX}.all.chain $target $query ${PREFIX}ChainPreNet.chain

echo "
chainNet ${PREFIX}ChainPreNet.chain $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin $net | netFilter /dev/stdout -chimpSyn > $net
"
chainNet ${PREFIX}ChainPreNet.chain $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin $net | netFilter /dev/stdout -chimpSyn > $net

echo "
netChainSubset $net ${PREFIX}ChainPreNet.chain ${PREFIX}.netted.chain
"
netChainSubset $net ${PREFIX}ChainPreNet.chain ${PREFIX}.netted.chain

echo "
netToAxt $net ${PREFIX}.all.chain $bitTarget $bitQuery $axt
"
netToAxt $net ${PREFIX}.all.chain $bitTarget $bitQuery $axt

echo "
axtDropOverlap $axt $target $query ${axt}.tmp; mv ${axt}.tmp $axt"
axtDropOverlap $axt $target $query ${axt}.tmp; mv ${axt}.tmp $axt

echo "
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $target $axt ${PREFIX}.sam; samtools sort -@ $SLURM_CPUS_ON_NODE ${PREFIX}.sam > ${PREFIX}.bam
"
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $target $axt ${PREFIX}.sam; samtools sort -@ $SLURM_CPUS_ON_NODE ${PREFIX}.sam > ${PREFIX}.bam

echo "
samtools index ${PREFIX}.bam"
samtools index ${PREFIX}.bam

rm ${PREFIX}.sam ${PREFIX}ChainPreNet.chain

echo "FINISHED"
