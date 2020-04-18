#!/bin/sh -e
#Either load modules or assign paths to GATK, samtools, and picard tools
module load GATK/4.1.3.0-gcb01

BAM=$1
PREFIX=$(basename $BAM .bam)

#Set your reference:
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
#if you want to name your files anything else set these variables

DIR=${PREFIX}"_genotypeVcfs"
mkdir -p $DIR
nodups=$DIR/${PREFIX}.nodups.bam
output=$DIR/${PREFIX}.gatk.valid.bam
gVcf=$DIR/${PREFIX}.wgs.g.vcf.gz

gatk --java-options "-Xmx4G" MarkDuplicates -I $BAM -O $nodups -M $DIR/${PREFIX}".dup.metrics" -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true
#adds Read group information to bam alignments
#this is how GATK differentiates between your cohort of samples

#required aruguments:
library="WGS"
platform="Illumina"
unit="HiSeqX"
sampleID=$PREFIX

gatk --java-options "-Xmx4G" AddOrReplaceReadGroups -I $nodups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $PREFIX | gatk --java-options "-Xmx4G" FixMateInformation -I /dev/stdin -O $output -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true
samtools index $output

#Now run GATK on processed BAM
gatk HaplotypeCaller --java-options "-Xmx16G" \
	-I $output -O $gVcf -R $REF \
	-ERC GVCF --native-pair-hmm-threads 8 \
rm $nodups ${nodups}.bai
