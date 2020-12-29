#!/bin/sh
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --time=0-02
#SBATCH --nodes=1
set -e

queryFa=gasAcu1-4.fa
targetFa=/data/lowelab/edotau/rabsTHREEspine/index/rabsTHREEspine.fa
CHAINQ=rabs3Bepa4Netted.chain

QNAME=$(basename $queryFa .fa)
TNAME=$(basename $targetFa .fa)

QCHROM=${QNAME}.sizes
TCHROM=${TNAME}.sizes

BITQUERY=${QNAME}.2bit
BITTARGET=${TNAME}.2bit

export PATH=/data/lowelab/edotau/kentUtils:/data/lowelab/edotau/bin/envs/htslib/bin:$PATH

echo "
chainStitchId $CHAINQ /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin ${QNAME}.${TNAME}.tBest.chain
"
chainStitchId $CHAINQ /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin ${QNAME}.${TNAME}.tBest.chain

echo "
chainPreNet ${QNAME}.${TNAME}.tBest.chain $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > ${QNAME}.${TNAME}.rBest.net.gz
"
chainPreNet ${QNAME}.${TNAME}.tBest.chain $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > ${QNAME}.${TNAME}.rBest.net.gz

echo "
netChainSubset ${QNAME}.${TNAME}.rBest.net.gz ${QNAME}.${TNAME}.tBest.chain /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > ${QNAME}.${TNAME}.rBest.netted.chain.gz
"
netChainSubset ${QNAME}.${TNAME}.rBest.net.gz ${QNAME}.${TNAME}.tBest.chain /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > ${QNAME}.${TNAME}.rBest.netted.chain.gz

echo "
chainSwap ${QNAME}.${TNAME}.rBest.netted.chain.gz /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}.${QNAME}.rBest.netted.chain.gz
"
chainSwap ${QNAME}.${TNAME}.rBest.netted.chain.gz /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}.${QNAME}.rBest.netted.chain.gz

echo "
chainPreNet ${TNAME}.${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout  | gzip -c > ${TNAME}.${QNAME}.rbest.net.gz
"
chainPreNet ${TNAME}.${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout  | gzip -c > ${TNAME}.${QNAME}.rbest.net.gz

rm ${QNAME}.${TNAME}.tBest.chain
md5sum *.rbest.*.gz > md5sum.rbest.txt

echo "
chainToPsl ${QNAME}.${TNAME}.rBest.netted.chain.gz $QCHROM $TCHROM $BITQUERY $BITTARGET ${QNAME}.${TNAME}.rBest.netted.chain.psl"
chainToPsl ${QNAME}.${TNAME}.rBest.netted.chain.gz $QCHROM $TCHROM $BITQUERY $BITTARGET ${QNAME}.${TNAME}.rBest.netted.chain.psl

echo "
chainToPsl ${TNAME}.${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM $BITTARGET $BITQUERY ${TNAME}.${QNAME}.rBest.netted.chain.psl"
chainToPsl ${TNAME}.${QNAME}.rBest.netted.chain.gz $TCHROM $QCHROM $BITTARGET $BITQUERY ${TNAME}.${QNAME}.rBest.netted.chain.psl

echo "
netToAxt ${TNAME}.${QNAME}.rbest.net.gz ${TNAME}.${QNAME}.rBest.netted.chain.gz $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}.${QNAME}.rBest.netted.chain.axt.gz
"
netToAxt ${TNAME}.${QNAME}.rbest.net.gz ${TNAME}.${QNAME}.rBest.netted.chain.gz $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > ${TNAME}.${QNAME}.rBest.netted.chain.axt.gz

echo "
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM ${TNAME}.${QNAME}.rBest.netted.chain.axt.gz ${TNAME}.${QNAME}.rBest.netted.chain.sam"
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM ${TNAME}.${QNAME}.rBest.netted.chain.axt.gz ${TNAME}.${QNAME}.rBest.netted.chain.sam

echo "
samtools sort -@ $SLURM_CPUS_ON_NODE ${TNAME}.${QNAME}.rBest.netted.chain.sam > ${TNAME}.${QNAME}.rBest.netted.chain.bam"
samtools sort -@ $SLURM_CPUS_ON_NODE ${TNAME}.${QNAME}.rBest.netted.chain.sam > ${TNAME}.${QNAME}.rBest.netted.chain.bam

echo "
samtools index ${TNAME}.${QNAME}.rBest.netted.chain.bam"
samtools index ${TNAME}.${QNAME}.rBest.netted.chain.bam

rm ${TNAME}.${QNAME}.rBest.netted.chain.sam
