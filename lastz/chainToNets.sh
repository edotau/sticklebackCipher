
#!/bin/sh
module add kentUtils/v302-gcb01
PREFIX=$1
contigs=$2
DIR=$3

queryPREFIX=$(basename $contigs .fa)
query=${queryPREFIX}.sizes
bitQuery=${queryPREFIX}.2bit
faToTwoBit $contigs $queryPREFIX
faSize -detailed $contigs > $query

target=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.sizes
bitTarget=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.2bit
faIdx=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta.fai

chainMergeSort $DIR*.chain | chainSplit chainMerge stdin -lump=50
cat chainMerge/*.chain > ${PREFIX}.all.chain
chainSort ${PREFIX}.all.chain ${PREFIX}.all.sorted.chain
rm -r chainMerge

chainPreNet ${PREFIX}.all.sorted.chain $target $query ${PREFIX}ChainPreNet.chain
net=${PREFIX}.net
chainNet ${PREFIX}ChainPreNet.chain $target $query /dev/stdout /dev/null | netSyntenic /dev/stdin $net
#${PREFIX}_syntenic_netted_chain.net
netFilter -minScore=500000 $net > ${PREFIX}Filtered.net

netChainSubset ${PREFIX}Filtered.net ${PREFIX}.all.sorted.chain ${PREFIX}Netted.chain
axt=${PREFIX}NettedChain.axt
netToAxt ${PREFIX}Filtered.net ${PREFIX}.all.sorted.chain $bitTarget $bitQuery $axt
axtToMaf $axt $target $query ${PREFIX}.chaining.maf

mafConvert=/data/lowelab/edotau/software/last-979/scripts/maf-convert
$mafConvert sam ${PREFIX}.chaining.maf > ${PREFIX}.tmp

module add samtools
samtools view -ht $faIdx ${PREFIX}.tmp | samtools view -h -b > ${PREFIX}.bam
samtools sort -T ${PREFIX}tmp.bam -o ${PREFIX}Chaining.bam ${PREFIX}.bam
samtools index ${PREFIX}Chaining.bam
