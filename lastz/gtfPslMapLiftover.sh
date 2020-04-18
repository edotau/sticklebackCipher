#!/bin/sh -e
module add gcc kentUtils
GTF=$1
chain=$2
query=$3
PREFIX=$(basename $GTF .gtf)
output=$(echo $query | cut -d '.' -f 1)
gp=${PREFIX}.gp
chromSizes=/data/lowelab/RefGenomes/gasAcu1/gasAcu1.sizes
/data/lowelab/cl454/bin/x86_64/gtfToGenePred -genePredExt -includeVersion $GTF $gp
genePredToPsl=/data/lowelab/cl454/bin/x86_64/genePredToPsl
psl=${PREFIX}.psl
$genePredToPsl $chromSizes $gp $psl
mappedPsl=${output}.psl
map=/data/lowelab/cl454/bin/x86_64/pslMap -chainMapFile -swapMap $psl $chain $mappedPsl
/data/lowelab/cl454/bin/x86_64/pslToBed $mappedPsl ${output}.bed
rm $gp $psl
bedToGenePred ${output}.bed ${output}.gp
genePredToGtf file ${output}.gtf
