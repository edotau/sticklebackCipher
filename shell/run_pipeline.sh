#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
set -e
module add GATK
one=$1
two=$2
BAM=$3
/data/lowelab/edotau/rabsTHREEspine/genePredMerge/rnaseq_bamfiles/mergeBam.sh $one $two $BAM
REF=/data/lowelab/edotau/rabsTHREEspine/index/refdata-rabsTHREEspine/fasta/genome.fa
PREFIX=$(basename $BAM .bam)
OUT=${PREFIX}.noSpanningSplicing.bam
/data/lowelab/ericAndCraig/gatk-4.1.8.1/gatk SplitNCigarReads --java-options "-Xmx12G" --input $BAM --output $OUT --reference $REF 2>${PREFIX}.SplitNCigarReads.log
module add samtools
samtools index $OUT
/data/lowelab/edotau/rabsTHREEspine/genePredMerge/rnaseq_bamfiles/rnaSeqCombine/bamToBigWig.sh $OUT
