#!/bin/sh
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --time=0-02
#SBATCH --nodes=1
set -e

export PATH=/data/lowelab/edotau/kentUtils/:/data/lowelab/edotau/bin/envs/htslib/bin/:$PATH

targetContigs=$1
queryContigs=$2
DIR=${3%/}

targetPREFIX=$(basename $targetContigs .fa)
queryPREFIX=$(basename $queryContigs .fa)

PREFIX=${targetPREFIX}_${queryPREFIX}

target=${targetPREFIX}.sizes
query=${queryPREFIX}.sizes

bitTarget=${targetPREFIX}.2bit
bitQuery=${queryPREFIX}.2bit

net=${PREFIX}.chain.net.gz
axt=${PREFIX}.netted.chain.axt

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
cat chainMerge/*.chain | grep -v '#' > ${PREFIX}.all.chain"
cat chainMerge/*.chain | grep -v '#' > ${PREFIX}.all.chain

echo "
chainSort ${PREFIX}.all.chain ${PREFIX}.sorted.chain"
chainSort ${PREFIX}.all.chain ${PREFIX}.sorted.chain
rm -r chainMerge

echo "
chainPreNet ${PREFIX}.sorted.chain $target $query /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > $net
"
chainPreNet ${PREFIX}.sorted.chain $target $query /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > $net

rm ${PREFIX}.all.chain
echo "
netChainSubset $net ${PREFIX}.sorted.chain ${PREFIX}.netted.chain
"
netChainSubset $net ${PREFIX}.sorted.chain ${PREFIX}.netted.chain

echo "
chainToAxt ${PREFIX}.netted.chain $bitTarget $bitQuery $axt"
chainToAxt ${PREFIX}.netted.chain $bitTarget $bitQuery $axt

echo "
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $target $axt ${PREFIX}.sam
"
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $target $axt ${PREFIX}.sam

echo "
samtools sort -@ $SLURM_CPUS_ON_NODE ${PREFIX}.sam > ${PREFIX}.bam"
samtools sort -@ $SLURM_CPUS_ON_NODE ${PREFIX}.sam > ${PREFIX}.bam

echo "
samtools index ${PREFIX}.bam"
samtools index ${PREFIX}.bam

netFileStats="netFileStats"
mkdir -p $netFileStats

echo "
netStats -gap=$netFileStats/${PREFIX}.gap.txt -fill=$netFileStats/${PREFIX}.fill.txt -top=$netFileStats/${PREFIX}.top.txt -syn=$netFileStats/${PREFIX}.syn.txt -nonSyn=$netFileStats/${PREFIX}.nonsyn.txt -syn=$netFileStats/${PREFIX}.syn.txt -inv=$netFileStats/${PREFIX}.inv.txt -dupe=$netFileStats/${PREFIX}.dupe.txt $netFileStats/${PREFIX}.summary.txt $net
"
netStats -gap=$netFileStats/${PREFIX}.gap.txt -fill=$netFileStats/${PREFIX}.fill.txt -top=$netFileStats/${PREFIX}.top.txt -syn=$netFileStats/${PREFIX}.syn.txt -nonSyn=$netFileStats/${PREFIX}.nonsyn.txt a-syn=$netFileStats/${PREFIX}.syn.txt -inv=$netFileStats/${PREFIX}.inv.txt -dupe=$netFileStats/${PREFIX}.dupe.txt $netFileStats/${PREFIX}.summary.txt $net

rm ${PREFIX}.sam

echo "FINISHED"
