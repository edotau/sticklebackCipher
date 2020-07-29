#!/bin/sh
#SBATCH --mem20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --job-name=SNPsToGenotypeVcf
#SBATCH --parsable
set -e
module load GATK/4.1.3.0-gcb01 samtools
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
snpDb=/data/lowelab/edotau/master/LITCxMATA/scratch/LITCxMATA.cohort.SNPs.recalVQSR.sitesonlyForBQSR.vcf.gz
indelDb=/data/lowelab/edotau/master/LITCxMATA/indels/scratch/LITCxMATA.cohort.GENOTYPED.INDEL.recalVQSR.sitesonlyForBQSR.vcf.gz
BAM=$1
DIR=$2
PREFIX=$(basename $BAM .bam)

#Set output variables:

DIR=recalGenotypeVcfs
recalBamDir=BAMSrecalb
mkdir -p $DIR
mkdir -p $recalBamDir
recalBam=$recalBamDir/${PREFIX}.recal.bam
recalTable=$recalBamDir/${PREFIX}.recal_data.table
genotypeVcf=$DIR/${PREFIX}.g.vcf.gz

gatk BaseRecalibrator -R $REF --known-sites $snpDb --known-sites  $indelDb -I $BAM -O $recalTable

gatk ApplyBQSR -R $REF -I $BAM --bqsr-recal-file $recalTable -O $recalBam

samtools index $recalBam
gatk HaplotypeCaller --input $recalBam --output $genotypeVcf --reference $REF -ERC GVCF
