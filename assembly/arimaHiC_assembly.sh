#!/bin/sh

REF=$1
READ1=$2
READ2=$3

PREFIX1=$(echo $READ1 | sed 's/.fastq.gz//')
PREFIX2=$(echo $READ2 | sed 's/.fastq.gz//')

bwa_filter=/data/lowelab/edotau/sticklebackCipher/assembly/bwa_filter.sh
ji1=$(sbatch $bwa_filter $READ1 $REF)
jid2=$(sbatch $bwa_filter $READ2 $REF)

mappingFilter=/data/lowelab/edotau/scripts/mappingQuality.sh
jid3=$(sbatch --dependency=afterany:$jid1:$jid2 $mappingFilter $REF ${PREFIX1}_filter.bam ${PREFIX2}_filter.bam)
