#!/bin/sh
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=0-06
#SBATCH --nodes=1
set -e
module add kentUtils/v302-gcb01
PREFIX=$1
targetContigs=$2
queryContigs=$3
DIR=${4%/}

targetPREFIX=$(basename $targetContigs .fa)
queryPREFIX=$(basename $queryContigs .fa)

target=${targetPREFIX}.sizes
query=${queryPREFIX}.sizes
faSize -detailed $targetContigs > $target
faSize -detailed $queryContigs > $query

bitTarget=${targetPREFIX}.2bit
bitQuery=${queryPREFIX}.2bit
faToTwoBit $targetContigs $bitTarget
faToTwoBit $queryContigs $bitQuery

chainMergeSort $DIR/*.chain | chainSplit chainMerge stdin -lump=50
cat chainMerge/*.chain > ${PREFIX}.all.chain
chainSort ${PREFIX}.all.chain ${PREFIX}.all.sorted.chain
rm -r chainMerge

chainPreNet ${PREFIX}.all.sorted.chain $target $query ${PREFIX}ChainPreNet.chain
net=${PREFIX}.net
chainNet ${PREFIX}ChainPreNet.chain $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin $net
#${PREFIX}_syntenic_netted_chain.net
netFilter -minScore=500000 $net > ${PREFIX}Filtered.net

netChainSubset ${PREFIX}Filtered.net ${PREFIX}.all.sorted.chain ${PREFIX}Netted.chain
axt=${PREFIX}.NettedChain.axt
netToAxt ${PREFIX}Filtered.net ${PREFIX}.all.sorted.chain $bitTarget $bitQuery $axt

axtToSam=/data/lowelab/edotau/RABS.HiCanuScaff10x/software/gonomics/cmd/axtSam/axtSam
$axtToSam -chrom $target $axt /dev/stdout | samtools sort -@ $SLURM_CPUS_ON_NODE > ${PREFIX}.NettedChain.bam
samtools index ${PREFIX}.NettedChain.bam
