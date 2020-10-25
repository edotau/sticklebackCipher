#!/bin/sh
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0-04
#SBATCH --nodes=1
set -e
module add bedtools2/2.25.0-fasrc01
module add coreutils/8.25-gcb01
module add kentUtils/v302-gcb01
BAM=$1
PREFIX=$(basename $BAM _final.bam)
bed=${PREFIX}.bed
chromSize=/data/lowelab/edotau/rabsTHREEspine/index/rabsTHREEspine.sizes
outBg=${PREFIX}.bg
BIGWIG=${PREFIX}.bw

echo "
bamToBed -i $BAM > $bed"
bamToBed -i $BAM > $bed

echo "
sort -k 1,1 $bed -T $PWD --parallel=$SLURM_CPUS_ON_NODE > ${PREFIX}.tmp && mv ${PREFIX}.tmp $bed"
sort -k 1,1 $bed -T $PWD --parallel=$SLURM_CPUS_ON_NODE > ${PREFIX}.tmp && mv ${PREFIX}.tmp $bed 

echo "
genomeCoverageBed -bg -i $bed -g $chromSize > $outBg"
genomeCoverageBed -bg -i $bed -g $chromSize > $outBg

echo "
bedGraphToBigWig $outBg $chromSize $BIGWIG"
bedGraphToBigWig $outBg $chromSize $BIGWIG

echo "
scp $BIGWIG eha17@trackhub.genome.duke.edu:/nfs/trackhub/lowelab/edotau/rabsTHREEspine/myHub/bams/RNASEQ"
scp $BIGWIG eha17@trackhub.genome.duke.edu:/nfs/trackhub/lowelab/edotau/rabsTHREEspine/myHub/bams/RNASEQ

echo "finished converting bed to bigwig! cleaning up tmp files:
rm $bed $outBg"
rm $bed $outBg
