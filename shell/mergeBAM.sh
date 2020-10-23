#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL --mail-user=eric.au@duke.edu
#SBATCH --job-name=RNASEQ.sh
set -e
module add samtools
OUT=$1
SORT_THREADS=$(($SLURM_CPUS_ON_NODE / 4))
THREADS=$(($SLURM_CPUS_ON_NODE - $SORT_THREADS))
echo "
samtools merge -@ $THREADS - *.bam | samtools sort -@ $SORT_THREADS -m 2G > $OUT"
samtools merge -@ $THREADS - *.bam | samtools sort -@ $SORT_THREADS -m 2G > $OUT

echo "
samtools index $OUT"
samtools index $OUT

echo "
$PWD/stringtie_exe.sh $OUT"
$PWD/stringtie_exe.sh $OUT
