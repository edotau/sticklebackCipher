#!/bin/sh
READ1=$1
READ2=$2
#module add fastqc
#fastqc -t 6 $READ1 $READ2
PREFIX=$(echo $READ1 | sed 's/_R1.fastq.gz//')
#echo "###### trimming reads, vis bbduk"
#/data/lowelab/software/bbmap/bbduk.sh \
#	in1=$READ1 \
#	in2=$READ2 \
#	out1=${PREFIX}_R1.trim.fastq.gz \
#	out2=${PREFIX}_R2.trim.fastq.gz \
#	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
#	ref=/data/lowelab/software/bbmap/resources/adapters.fa

module load cutadapt/1.8.3-gcb01
cutadapt -a CTGTCTCTTATACACATCT -o ${PREFIX}_R1.trim.fastq.gz $READ1
cutadapt -a CTGTCTCTTATACACATCT -o ${PREFIX}_R2.trim.fastq.gz $READ2



module load bwa samtools
module load bedtools2/2.25.0-fasrc01
module load kentUtils/v302-gcb01

mREF=/data/lowelab/RefGenomes/rabs/gasAcu2RABS.fasta
mSize=/data/lowelab/RefGenomes/rabs/gasAcu2RABS.sizes
#fRef=/data/lowelab/RefGenomes/gasAcu1/gasAcu1.fa
#fSize=/data/lowelab/RefGenomes/gasAcu1/gasAcu1.sizes

############PEAK CALLER#############
callPeaks=/data/lowelab/software/Genrich-master/Genrich

marineDIR=gasAcu2RABS_12_09_19
#freshDIR=fresh

mkdir -p $marineDIR
#mkdir -p $freshDIR
echo "Starting ATACSEQ Analysis for Marine genome"
bwa mem -t 12 $mREF ${PREFIX}_R1.trim.fastq.gz ${PREFIX}_R2.trim.fastq.gz  > ${marineDIR}/${PREFIX}.sam
samtools view -b -h -q 30 ${marineDIR}/${PREFIX}.sam > ${marineDIR}/${PREFIX}.bam
samtools sort -n -o ${marineDIR}/${PREFIX}.name.bam ${marineDIR}/${PREFIX}.bam
samtools fixmate -m ${marineDIR}/${PREFIX}.name.bam ${marineDIR}/${PREFIX}.fixmate.bam
samtools sort -T ${marineDIR}/${PREFIX}_tmp -o ${marineDIR}/${PREFIX}.sorted.bam ${marineDIR}/${PREFIX}.fixmate.bam
samtools markdup -r ${marineDIR}/${PREFIX}.sorted.bam ${marineDIR}/${PREFIX}.sorted.nodups.bam
samtools index ${marineDIR}/${PREFIX}.sorted.nodups.bam
samtools stats ${marineDIR}/${PREFIX}.sorted.nodups.bam > ${marineDIR}/${PREFIX}.stats
echo "Finished mapping and marking dups"
rm ${marineDIR}/${PREFIX}.sam ${marineDIR}/${PREFIX}.name.bam ${marineDIR}/${PREFIX}.fixmate.bam ${marineDIR}/${PREFIX}.sorted.bam
#call peaks
marinePeaks=${marineDIR}/peakFiles
mkdir -p $marinePeaks

samtools sort -n ${marineDIR}/${PREFIX}.sorted.nodups.bam -o ${marineDIR}/${PREFIX}.4Peaks.bam
$callPeaks -t ${marineDIR}/${PREFIX}.4Peaks.bam -o ${marinePeaks}/${PREFIX}.bed -j  -y  -r  -e chrM  -v -a 25 -q 0.05
sort -k1,1 -k2,2n ${marinePeaks}/${PREFIX}.bed > ${marinePeaks}/${PREFIX}.sorted.bed

echo "Finished Calling Peaks"
genomeCoverageBed -bg -i ${marinePeaks}/${PREFIX}.sorted.bed -g $mSize > ${marinePeaks}/${PREFIX}.bg
bedGraphToBigWig ${marinePeaks}/${PREFIX}.bg $mSize ${marinePeaks}/${PREFIX}.bw

#echo "Starting ATACSEQ Analysis for Freshwater genome"
#bwa mem -t 12 $fREF ${PREFIX}_R1.trim.fastq.gz ${PREFIX}_R2.trim.fastq.gz  > ${freshDIR}/${PREFIX}.sam
#samtools view -b -h -q 30 ${freshDIR}/${PREFIX}.sam > ${freshDIR}/${PREFIX}.bam
#samtools sort -n -o ${freshDIR}/${PREFIX}.name.bam ${freshDIR}/${PREFIX}.bam
#samtools fixmate -m ${freshDIR}/${PREFIX}.name.bam ${freshDIR}/${PREFIX}.fixmate.bam
#samtools sort -T ${freshDIR}/${PREFIX}_tmp -o ${freshDIR}/${PREFIX}.sorted.bam ${freshDIR}/${PREFIX}.fixmate.bam
#samtools markdup -r ${freshDIR}/${PREFIX}.sorted.bam ${freshDIR}/${PREFIX}.sorted.nodups.bam
#samtools index ${freshDIR}/${PREFIX}.sorted.nodups.bam
#samtools stats ${freshDIR}/${PREFIX}.sorted.nodups.bam > ${freshDIR}/${PREFIX}.stats
#echo "Finished mapping and marking dups"
#rm ${freshDIR}/${PREFIX}.sam ${freshDIR}/${PREFIX}.name.bam ${freshDIR}/${PREFIX}.fixmate.bam ${freshDIR}/${PREFIX}.sorted.bam
#call peaks
#freshPeaks=${freshDIR}/peakFiles
#mkdir -p $freshPeaks

#samtools sort -n ${freshDIR}/${PREFIX}.sorted.nodups.bam -o ${freshDIR}/${PREFIX}.sorted
#$callPeaks -t ${freshDIR}/${PREFIX}.sorted -o ${freshPeaks}/${PREFIX}.bed -j  -y  -r  -e chrM  -v -a 25 -q 0.01
#sort -k1,1 -k2,2n ${freshPeaks}/${PREFIX}.bed > ${freshPeaks}/${PREFIX}.sorted.bed

#echo "Finished Calling Peaks"
#genomeCoverageBed -bg -i ${freshPeaks}/${PREFIX}.sorted.bed -g $fSize > ${freshPeaks}/${PREFIX}.bg
#bedGraphToBigWig ${freshPeaks}/${PREFIX}.bg $fSize ${freshPeaks}/${PREFIX}.bw

rm ${PREFIX}_R1.trim.fastq.gz
rm ${PREFIX}_R2.trim.fastq.gz

