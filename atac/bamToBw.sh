#!/bin/sh
module load bwa samtools
module load bedtools2/2.25.0-fasrc01
module load kentUtils/v302-gcb01
module add bedtools2/2.25.0-fasrc01
SIZE=/data/lowelab/RefGenomes/gasAcu1/gasAcu1.sizes
REF=/data/lowelab/RefGenomes/gasAcu1/gasAcu1.fa
BAM=$1
PREFIX=$(echo $BAM | sed 's/.bam//')
bed=${PREFIX}.bed
bamToBed -i $BAM > $bed
genomeCoverageBed -bg -i $bed -g $SIZE > ${PREFIX}.bg
bedGraphToBigWig ${PREFIX}.bg $SIZE ${PREFIX}.bw
rm ${PREFIX}.bg
