#!/bin/sh -e
GTF=$1
chromSizes=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.sizes
PREFIX=$(basename $GTF .gtf)
module load gcc
genePred=${PREFIX}.gp
/data/lowelab/cl454/bin/x86_64/gtfToGenePred -genePredExt -includeVersion $GTF $genePred
/data/lowelab/cl454/bin/x86_64/genePredToBigGenePred -known $genePred stdout | sort -k1,1 -k2,2n > ${PREFIX}.bgp
/data/lowelab/cl454/bin/x86_64/bedToBigBed -type=bed12+8 -tab -as=/data/lowelab/edotau/scripts/bigGenePred.as ${PREFIX}.bgp $chromSizes ${PREFIX}.bb
rm $genePred
echo "FINISHED!"
