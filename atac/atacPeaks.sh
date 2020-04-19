#!/bin/sh
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --time=1-00:00:00
BAM=$1
PREFIX=$(basename $BAM .bam)
qNameSort=${PREFIX}.name.sorted.bam
samtools sort -@ 8 -n -o $qNameSort $BAM
q=(0.05 0.01)
a=(200 300)
picard=/data/lowelab/edotau/software/picard.jar
module add samtools java/1.8.0_45-fasrc01
java -jar $picard FixMateInformation I=$qNameSort ADD_MATE_CIGAR=true IGNORE_MISSING_MATES=true
peakCaller=/data/lowelab/software/Genrich-master/Genrich
$peakCaller -t $qNameSort -o ${PREFIX}.tmp -j -y -r -v -D
sort -k1,1 -k2,2n ${PREFIX}.tmp > ${PREFIX}.bed
