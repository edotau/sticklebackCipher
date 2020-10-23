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

module add StringTie/2.1.1-gcb01
#GTF=/data/lowelab/edotau/rabsTHREEspine/index/bepa_v1_ensembl_pslMappedRabs.gtf
#BAM=$1
#PREFIX=$(basename $BAM .bam)
#stringtie $BAM -p $SLURM_CPUS_ON_NODE -o ${PREFIX}.gtf -l RNASEQ
