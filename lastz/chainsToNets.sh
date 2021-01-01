#!/bin/sh
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --time=0-02
#SBATCH --nodes=1
set -e

export PATH=/data/lowelab/edotau/kentUtils/:/data/lowelab/edotau/bin/envs/htslib/bin/:$PATH

TARGET=$1
QUERY=$2
DIR=${3%/}

TNAME=$(basename $TARGET .fa)
QNAME=$(basename $QUERY .fa)

PREFIX=${TNAME}_${QNAME}

TCHROM=${TNAME}.sizes
QCHROM=${QNAME}.sizes

BITTARGET=${TNAME}.2bit
BITQUERY=${QNAME}.2bit

NET=${PREFIX}.chain.net
AXT=${PREFIX}.netted.chain.axt

if [[ "$#" -lt 3 ]]
then
	echo "Usage: ./chainsToNets.sh target.fa query.fa /path/to/chains"
	exit 0
else
	if ! [ -e "$BITTARGET" ] ; then
		echo "faToTwoBit $TARGET $BITTARGET
		"
        	faToTwoBit $TARGET $BITTARGET
	fi

	if ! [ -e "$BITQUERY" ] ; then
		echo "faToTwoBit $QUERY $BITQUERY
		"
        	faToTwoBit $QUERY $BITQUERY
	fi

	if ! [ -e "$TCHROM" ] ; then
		echo "faSize -detailed $TARGET > $TCHROM
		"
        	faSize -detailed $TARGET > $TCHROM
	fi

	if ! [ -e "$QCHROM" ] ; then
		echo "faSize -detailed $QUERY > $QCHROM
		"
        	faSize -detailed $QUERY > $QCHROM
	fi
fi

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
chainPreNet ${PREFIX}.sorted.chain $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin $NET
"
chainPreNet ${PREFIX}.sorted.chain $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin $NET

rm ${PREFIX}.all.chain
echo "
netChainSubset $NET ${PREFIX}.sorted.chain ${PREFIX}.netted.chain
"
netChainSubset $NET ${PREFIX}.sorted.chain ${PREFIX}.netted.chain

echo "
chainToAxt ${PREFIX}.netted.chain $BITTARGET $BITQUERY $AXT"
chainToAxt ${PREFIX}.netted.chain $BITTARGET $BITQUERY $AXT

echo "
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM $AXT ${PREFIX}.sam
"
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM $AXT ${PREFIX}.sam

echo "
samtools sort -@ $SLURM_CPUS_ON_NODE ${PREFIX}.sam > ${PREFIX}.bam"
samtools sort -@ $SLURM_CPUS_ON_NODE ${PREFIX}.sam > ${PREFIX}.bam

echo "
samtools index ${PREFIX}.bam"
samtools index ${PREFIX}.bam

netFileStats="netFileStats"
mkdir -p $netFileStats

echo "
netStats -gap=$netFileStats/${PREFIX}.gap.txt -fill=$netFileStats/${PREFIX}.fill.txt -top=$netFileStats/${PREFIX}.top.txt -syn=$netFileStats/${PREFIX}.syn.txt -nonSyn=$netFileStats/${PREFIX}.nonsyn.txt -syn=$netFileStats/${PREFIX}.syn.txt -inv=$netFileStats/${PREFIX}.inv.txt -dupe=$netFileStats/${PREFIX}.dupe.txt $netFileStats/${PREFIX}.summary.txt $NET
"
netStats -gap=$netFileStats/${PREFIX}.gap.txt -fill=$netFileStats/${PREFIX}.fill.txt -top=$netFileStats/${PREFIX}.top.txt -syn=$netFileStats/${PREFIX}.syn.txt -nonSyn=$netFileStats/${PREFIX}.nonsyn.txt a-syn=$netFileStats/${PREFIX}.syn.txt -inv=$netFileStats/${PREFIX}.inv.txt -dupe=$netFileStats/${PREFIX}.dupe.txt $netFileStats/${PREFIX}.summary.txt $NET

rm ${PREFIX}.sam

echo "FINISHED"

