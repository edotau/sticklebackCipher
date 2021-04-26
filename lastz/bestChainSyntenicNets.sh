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
QRY_TGT=${QNAME}.${TNAME}.tBest.chain

echo "
chainStitchId $TCHAIN /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin $QRY_TGT
"
chainStitchId $TCHAIN /dev/stdout | chainSwap /dev/stdin /dev/stdout | chainSort /dev/stdin $QRY_TGT

# Net those on QUERY-REFERENCED to get reciprocal best net ####

QRY_TGT_NET=${QNAME}.${TNAME}.rBest.net.gz
echo "
chainPreNet $QRY_TGT $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > $QRY_TGT_NET
"
chainPreNet $QRY_TGT $QCHROM $TCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $QCHROM $TCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout | gzip -c > $QRY_TGT_NET

# Extract QUERY-REFERENCED reciprocal best chain ####

QRY_TGT_NETCHAIN=${QNAME}.${TNAME}.rBest.netted.chain.gz
echo "
netChainSubset $QRY_TGT_NET $QRY_TGT /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > $QRY_TGT_NETCHAIN
"
netChainSubset $QRY_TGT_NET $QRY_TGT /dev/stdout | chainStitchId /dev/stdin /dev/stdout | gzip -c > $QRY_TGT_NETCHAIN

# Swap to get TARGET-REFERENCED reciprocal best chain ####
TGT_QRY_SWAP_CHAIN=${TNAME}.${QNAME}.rBest.netted.chain.gz
echo "
chainSwap $QRY_TGT_NETCHAIN /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_CHAIN
"
chainSwap $QRY_TGT_NETCHAIN /dev/stdout | chainSort /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_CHAIN

# Net those on TARGET to get TARGET-REFERENCE reciprocal best net ###

TGT_QRY_SWAP_NET=${TNAME}.${QNAME}.rbest.net.gz
echo "
chainPreNet $TGT_QRY_SWAP_CHAIN $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout  | gzip -c > $TGT_QRY_SWAP_NET
"
chainPreNet $TGT_QRY_SWAP_CHAIN $TCHROM $QCHROM /dev/stdout | chainNet -minSpace=1 -minScore=0 /dev/stdin $TCHROM $QCHROM /dev/stdout /dev/null | netSyntenic /dev/stdin /dev/stdout  | gzip -c > $TGT_QRY_SWAP_NET

rm $QRY_TGT
md5sum *.rbest.*.gz > md5sum.rbest.txt

TGT_QRY_SWAP_NET_AXT=${TNAME}.${QNAME}.rBest.netted.chain.axt.gz
echo "
netToAxt $TGT_QRY_SWAP_NET ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_NET_AXT
"
netToAxt $TGT_QRY_SWAP_NET ${TNAME}THREE${QNAME}.rBest.netted.chain.gz $BITTARGET $BITQUERY /dev/stdout | axtSort /dev/stdin /dev/stdout | gzip -c > $TGT_QRY_SWAP_NET_AXT

AXT_SAM=${TNAME}.${QNAME}.rBest.netted.chain.sam
AXT_BAM=${TNAME}.${QNAME}.rBest.netted.chain.bam

echo "
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM $TGT_QRY_SWAP_NET_AXT $AXT_SAM
"
/data/lowelab/edotau/golang/src/github.com/vertgenlab/gonomics/cmd/axtSam/axtSam -chrom $TCHROM $TGT_QRY_SWAP_NET_AXT $AXT_SAM
echo "
samtools sort -@ $SLURM_CPUS_ON_NODE $AXT_SAM > $AXT_BAM"
samtools sort -@ $SLURM_CPUS_ON_NODE $AXT_SAM > $AXT_BAM"


echo "
samtools index $AXT_BAM
"
samtools index $AXT_BAM

rm $AXT_SAM 

