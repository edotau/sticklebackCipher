#!/bin/sh
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --job-name=SNPs
module load GATK/4.1.3.0-gcb01
BAM=$1
PREFIX=$(basename $BAM .gatk.valid.recal.bam)
DIR=GVCFs
mkdir -p $DIR
GVCF=$DIR/${PREFIX}.g.vcf.gz
gatk HaplotypeCaller --java-options "-Xmx50G" --input $BAM --output $GVCF --reference /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation --native-pair-hmm-threads 12
