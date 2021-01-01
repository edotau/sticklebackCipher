#!/bin/sh
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --time=1-0
#SBATCH --nodes=1
set -e

# export kentUils and samtools htslib to path or set the variables yourself
export PATH=/data/lowelab/edotau/kentUtils:/data/lowelab/edotau/bin/envs/htslib/bin:$PATH
# target fasta file
TARGET=$1

# query fasta file
QUERY=$2

# CHAINS target referenced
TCHAIN=$3

TNAME=$(basename $TARGET .fa)
QNAME=$(basename $QUERY .fa)

BITTARGET=${TNAME}.2bit
BITQUERY=${QNAME}.2bit

TCHROM=${TNAME}.sizes
QCHROM=${QNAME}.sizes

if [[ "$#" -lt 3 ]]
then
	echo "Usage: ./chainBestSynNets.sh target.fa query.fa target.chain"
	exit 1
else
	if ! [ -e "${TNAME}.2bit" ] ; then
        	faToTwoBit $TARGET $BITTARGET
	fi

	if ! [ -e "${QNAME}.2bit" ] ; then
        	faToTwoBit $QUERY $BITQUERY
	fi

	if ! [ -e "${TNAME}.sizes" ] ; then
        	faSize -detailed $TARGET > $TCHROM
	fi

	if ! [ -e "${QNAME}.2bit" ] ; then
        	faSize -detailed $QUERY > $QCHROM
	fi
fi

###### Swap TARGET best chains to be QUERY-REFERENCED ####
echo "
chainStitchId $TCHAIN /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin ${QNAME}THREE${TNAME}.tBest.chain
"
chainStitchId $TCHAIN /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin ${QNAME}THREE${TNAME}.tBest.chain

# Net those on QUERY-REFERENCED to get reciprocal best net ####
echo "
chainPreNet ${QNAME}THREE${TNAME}.tBest.chain $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > ${QNAME}THREE${TNAME}.rBest.net.gz
"
chainPreNet ${QNAME}THREE${TNAME}.tBest.chain $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > ${QNAME}THREE${TNAME}.rBest.net.gz

# Extract QUERY-REFERENCED reciprocal best chain ####
echo "
netChainSubset ${QNAME}THREE${TNAME}.rBest.net.gz ${QNAME}THREE${TNAME}.tBest.chain /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > ${QNAME}THREE${TNAME}.rBest.netted.chain.gz
"
netChainSubset ${QNAME}THREE${TNAME}.rBest.net.gz ${QNAME}THREE${TNAME}.tBest.chain /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > ${QNAME}THREE${TNAME}.rBest.netted.chain.gz

# Swap to get TARGET-REFERENCED reciprocal best chain ####
echo "
chainSwap ${QNAME}THREE${TNAME}.rBest.netted.chain.gz /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}THREE${QNAME}.rBest.netted.chain.gz
"
chainSwap ${QNAME}THREE${TNAME}.rBest.netted.chain.gz /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}THREE${QNAME}.rBest.netted.chain.gz

# Net those on TARGET to get TARGET-REFERENCE reciprocal best net ###
echo "
chainPreNet ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout  | gzip -c > ${TNAME}THREE${QNAME}.rbest.net.gz
"
chainPreNet ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout  | gzip -c > ${TNAME}THREE${QNAME}.rbest.net.gz

rm ${QNAME}THREE${TNAME}.tBest.chain
md5sum *.rbest.*.gz > md5sum.rbest.txt

echo "
netToBed -maxGap=1 ${QNAME}THREE${TNAME}.rBest.net.gz ${QNAME}THREE${TNAME}.rBest.net.bed"
netToBed -maxGap=1 ${QNAME}THREE${TNAME}.rBest.net.gz ${QNAME}THREE${TNAME}.rBest.net.bed

echo "
netToBed -maxGap=1 ${TNAME}THREE${QNAME}.rbest.net.gz ${TNAME}THREE${QNAME}.rbest.net.bed"
netToBed -maxGap=1 ${TNAME}THREE${QNAME}.rbest.net.gz ${TNAME}THREE${QNAME}.rbest.net.bed

echo "
chainToPsl ${QNAME}THREE${TNAME}.rBest.netted.chain.gz $QCHROM $TCHROM $BITQUERY $BITTARGET ${QNAME}THREE${TNAME}.rBest.netted.chain.psl"
chainToPsl ${QNAME}THREE${TNAME}.rBest.netted.chain.gz $QCHROM $TCHROM $BITQUERY $BITTARGET ${QNAME}THREE${TNAME}.rBest.netted.chain.psl

echo "
chainToPsl ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM $BITTARGET $BITQUERY ${TNAME}THREE${QNAME}.rBest.netted.chain.psl"
chainToPsl ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM $BITTARGET $BITQUERY ${TNAME}THREE${QNAME}.rBest.netted.chain.psl

echo "
netToAxt ${TNAME}.${QNAME}.rbest.net.gz ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}THREE${QNAME}.rBest.netted.chain.axt.gz
"
netToAxt ${TNAME}THREE${QNAME}.rbest.net.gz ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}THREE${QNAME}.rBest.netted.chain.axt.gz

echo "
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM ${TNAME}THREE${QNAME}.rBest.netted.chain.axt.gz ${TNAME}THREE${QNAME}.rBest.netted.chain.sam"
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM ${TNAME}THREE${QNAME}.rBest.netted.chain.axt.gz ${TNAME}THREE${QNAME}.rBest.netted.chain.sam

echo "
samtools sort -@ $SLURM_CPUS_ON_NODE ${TNAME}THREE${QNAME}.rBest.netted.chain.sam > ${TNAME}THREE${QNAME}.rBest.netted.chain.bam"
samtools sort -@ $SLURM_CPUS_ON_NODE ${TNAME}THREE${QNAME}.rBest.netted.chain.sam > ${TNAME}THREE${QNAME}.rBest.netted.chain.bam

echo "
samtools index ${TNAME}THREE${QNAME}.rBest.netted.chain.bam"
samtools index ${TNAME}THREE${QNAME}.rBest.netted.chain.bam

rm ${TNAME}THREE${QNAME}.rBest.netted.chain.sam

