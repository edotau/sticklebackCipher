#!/bin/sh
module load samtools bwa perl/5.10.1-fasrc04 java/1.8.0_45-fasrc01

bwa_filter=/data/lowelab/edotau/scripts/bwa_filter.sh
REF=$1
READ1=$2
READ2=$3

PREFIX1=$(echo $READ1 | sed 's/.fastq.gz//')
PREFIX2=$(echo $READ2 | sed 's/.fastq.gz//')

jid1=$(sbatch \
	--mem=64G \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=8 \
	--mail-type=ALL \
	--mail-user=eric.au@duke.edu \
	--job-name=hic_bwa_mapping01 \
	$bwa_filter $READ1 $REF)

jid2=$(sbatch \
	--mem=64G \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=8 \
	--mail-type=ALL \
	--mail-user=eric.au@duke.edu \
	--job-name=hic_bwa_mapping02 \
	$bwa_filter $READ2 $REF)

mappingFilter=/data/lowelab/edotau/scripts/mappingQuality.sh

jid3=$(sbatch --dependency=afterany:$jid1:$jid2 \
	--mem=64G \
	--nodes=1 \
	--ntasks=2 \
	--cpus-per-task=8 \
	--mail-type=ALL \
	--mail-user=eric.au@duke.edu \
	$mappingFilter $REF ${PREFIX1}_filter.bam ${PREFIX2}_filter.bam)
