#!/bin/sh
BAM=$1
PREFIX=$(echo $1 | sed 's/.sorted.nodups.bam.sorted//')
q=(0.05 0.01)
a=(200, 300)
Genrich=/data/lowelab/software/Genrich-master/Genrich
chromsizes=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.sizes
module load kentUtils
mkdir -p bedFiles
mkdir -p pileup
samtools sort -n ${BAM} -o ${BAM}.sorted
for i in "${q[@]}"
do
        for j in "${a[@]}"
        do
               
                $Genrich -t ${BAM}.sorted -o ${PREFIX}_${i}_${j}.tmp -j -q $i -a $j -k pileup/${PREFIX}_${i}_${j}.pileups.bg
                sort -k1,1 -k2,2n ${PREFIX}_${i}_${j}.tmp > bedFiles/${PREFIX}_${i}_${j}.bed
                echo $i $j
                rm ${PREFIX}_${i}_${j}.tmp
        done
done
exit 0

