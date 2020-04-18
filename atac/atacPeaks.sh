#!/bin/sh
BAM=$1
PREFIX=$(basename $BAM .sorted.bam.sorted)Pvalue
q=(0.05 0.01)
a=(200 300)
picard=/data/lowelab/edotau/software/picard.jar
module add samtools java/1.8.0_45-fasrc01
java -jar $picard FixMateInformation I=$BAM ADD_MATE_CIGAR=true IGNORE_MISSING_MATES=true ASSUME_SORTED=true
Genrich=/data/lowelab/software/Genrich-master/Genrich
$Genrich -t $BAM -o ${PREFIX}.tmp -j -y -r -e chrM -v
sort -k1,1 -k2,2n ${PREFIX}.tmp > ${PREFIX}.bed
