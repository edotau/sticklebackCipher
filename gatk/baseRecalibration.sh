#!/bin/sh
BAM=$1
PREFIX=$(basename $BAM .bam)
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --job-name=SNPsToGenotypeVcf
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=eric.au@duke.edu
module load GATK/4.1.3.0-gcb01 samtools
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
snpDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDbLITCxMATAGenotypeGVCFs.Filtered.SNP.vcf.gz
indelDb=/data/lowelab/edotau/toGasAcu2RABS/gVcf_LITXxMATA/filterDbLITCxMATAGenotypeGVCFs.Filtered.INDEL.vcf.gz
BAM=$1
PREFIX=$(basename $BAM .bam)

#Set output variables:
DIR=${PREFIX}_recalGenotypeVcfs
mkdir -p $DIR
recalBam=$DIR/${PREFIX}.recal.bam
recalTable=$DIR/${PREFIX}.recal_data.table
genotypeVcf=$DIR/${PREFIX}.g.vcf.gz

#NOTE: second known site argument is to mask certain varience
gatk --java-options "-Xmx32g" BaseRecalibrator -R $REF \
	--known-sites $snpDb \
       	--known-sites  $indelDb \
	-I $BAM -O $recalTable

gatk --java-options "-Xmx32g" ApplyBQSR -R $REF \
	-I $BAM --bqsr-recal-file $recalTable \
	-O $recalBam

samtools index $recalBam
gatk HaplotypeCaller --java-options "-Xmx32G" \
	--input $recalBam --output $genotypeVcf \
	--reference $REF \
	-ERC GVCF --native-pair-hmm-threads 12
