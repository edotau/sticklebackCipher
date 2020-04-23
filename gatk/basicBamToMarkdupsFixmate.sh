#!/bin/sh -e
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=toGenotypeVcf

#Either load modules or assign paths to GATK, samtools, and picard tools
module load GATK/4.1.3.0-gcb01
#Set your reference:
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta

BAM=$1
PREFIX=$(basename $BAM .bam)
#if you want to name your files anything else set these variables

DIR=$2
mkdir -p $DIR
markedDups=$DIR/${PREFIX}.markedDups.bam
output=$DIR/${PREFIX}.gatk.valid.bam
#final output
gVcf=${PREFIX}.g.vcf.gz

gatk --java-options "-Xmx16G" MarkDuplicates -I $BAM -O $markedDups -M $DIR/${PREFIX}".dup.metrics" -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true
#adds Read group information to bam alignments
#this is how GATK differentiates between your cohort of samples
#required aruguments:
library="WGS"
platform="Illumina"
unit="HiSeqX"
sampleID=$PREFIX

samtools view -b -@ 8 -f 1 $markedDups | gatk --java-options "-Xmx16G" AddOrReplaceReadGroups -I /dev/stdin -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $PREFIX | gatk --java-options "-Xmx16G" FixMateInformation -I /dev/stdin -O $output -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR $DIR/"fixmate.tmp"
samtools index $output
echo "Ready to run GATK on processed $PREIX BAM"
