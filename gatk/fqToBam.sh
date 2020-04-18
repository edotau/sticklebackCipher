#!/bin/sh

REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
module load bwa samtools/1.9-gcb01 GATK/4.1.3.0-gcb01
#used by trim_galore
module load cutadapt/2.3-gcb01 python/3.7.4-gcb01 pigz/2.3.4-gcb01

READ1=$1
READ2=$2

PREFIX=$(basename $READ1 _1.fq.gz | cut -d '_' -f 1)
DIR=${PREFIX}"to"$(basename $REF | cut -d '.' -f 1)
mkdir -p $DIR

trim_galore=/data/lowelab/edotau/software/TrimGalore-0.6.5/trim_galore
$trim_galore -o $DIR --basename $PREFIX --trim-n --max_n 0 --cores 4 --paired $READ1 $READ2
trim1=$DIR/${PREFIX}_val_1.fq.gz
trim2=$DIR/${PREFIX}_val_2.fq.gz

input1=$DIR/${PREFIX}_R1.fastq.gz
input2=$DIR/${PREFIX}_R2.fastq.gz

#Second Trimming of adaptors:
/data/lowelab/software/bbmap/bbduk.sh \
	in1=$trim1 \
	in2=$trim2 \
	out1=$input1 \
	out2=$input2 \
	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
	ref=/data/lowelab/software/bbmap/resources/adapters.fa
#MAPPING

echo -e "["$(date)"]\tAligning and sorting..."
BAM=$DIR/${PREFIX}.bam
bwa mem -t 8 $REF $1 $2 | samtools view -hb /dev/stdin | samtools sort -@ 4 -T $DIR/${PREFIX}.tmp -o $BAM
samtools index $BAM
