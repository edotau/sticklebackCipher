#!/bin/sh -e
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --time=1-00:00:00
module load bwa samtools
module load bedtools2/2.25.0-fasrc01
module load kentUtils/v302-gcb01
BAM=$1
REF=$2
mkdir -p $DIR
DIR=${PREFIX}toPeakCalls
PREFIX=$(echo $BAM | sed 's/.bam//')
bed=${PREFIX}.bed
bamToBed -i $BAM > $bed
#chromSize=$(echo $REF | cut -d '.' -f 1).sizes
faSize -detailed $REF | genomeCoverageBed -bg -i $bed -g /dev/stdin > ${PREFIX}.bg
bedGraphToBigWig ${PREFIX}.bg $SIZE ${PREFIX}.bw
rm ${PREFIX}.bg $bed
