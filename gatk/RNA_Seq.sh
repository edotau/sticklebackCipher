#!/bin/sh -e
module load GATK/4.1.0.0-gcb01 java/1.8.0_45-fasrc01 samtools/1.9-gcb01 STAR/2.7.2b-gcb01
picard=/data/lowelab/edotau/software/picard.jar

#used by trim_galore
module load cutadapt/2.3-gcb01 python/3.7.4-gcb01 pigz/2.3.4-gcb01

star_ref=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/STARidx

REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
input1=$1
input2=$2

PREFIX=$(echo $input1 | sed 's/w1_L8_R1_val_1.fq.gz//')toRABS
DIR=$PREFIX"rnaseq"
mkdir -p $DIR

trim_galore=/data/lowelab/edotau/software/TrimGalore-0.6.5/trim_galore
$trim_galore -o $DIR --basename $PREFIX --trim-n --max_n 0 --cores 4 --paired $input1 $input2

trim1=$DIR/${PREFIX}_val_1.fq.gz
trim2=$DIR/${PREFIX}_val_2.fq.gz

READ1=$DIR/${PREFIX}_R1.fastq.gz
READ2=$DIR/${PREFIX}_R2.fastq.gz

#Preprocessing
#Trimming adaptors
/data/lowelab/software/bbmap/bbduk.sh \
	in1=$trim1 \
	in2=$trim2 \
	out1=$READ1 \
	out2=$READ2 \
	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
	ref=/data/lowelab/software/bbmap/resources/adapters.fa

#MAPPING

echo -e "["$(date)"]\tAligning.."
STAR --outFileNamePrefix $DIR/$PREFIX --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMattrRGline ID:$PREFIX CN:Gonomics LB:PairedEnd PL:Illumina PU:Unknown SM:$PREFIX --genomeDir $star_ref --runThreadN 8 --readFilesCommand zcat --readFilesIn $READ1 $READ2 --twopassMode Basic

echo -e "["$(date)"]\tSorting.."
samtools sort -o $DIR/$PREFIX"_sorted.bam" $DIR/$PREFIX"Aligned.out.bam"
rm $DIR/$PREFIX"Aligned.out.bam"

echo -e "["$(date)"]\tIndexing.."
samtools index $DIR/$PREFIX"_sorted.bam"

echo -e "["$(date)"]\tMarking duplicates.."
java -Xmx8G -jar $picard MarkDuplicates I=$DIR/$PREFIX"_sorted.bam" O=$DIR/$PREFIX"_dupMarked.bam" M=$DIR/$PREFIX"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$DIR/$PREFIX.MarkDuplicates.log
rm $DIR/$PREFIX"_sorted.bam"
rm $DIR/$PREFIX"_sorted.bam.bai"

#SplitNCigarReads
echo -e "["$(date)"]\tSpliting reads.."
gatk SplitNCigarReads --java-options "-Xmx8g" --input $DIR/$PREFIX"_dupMarked.bam" --output $DIR/$PREFIX"_split.bam" --reference $REF 2>$DIR/$PREFIX.SplitNCigarReads.log
samtools index $DIR/$PREFIX"_split.bam"


gatk HaplotypeCaller --java-options "-Xmx8g" --input $DIR/$PREFIX"_split.bam" --dont-use-soft-clipped-bases true --output $DIR/${PREFIX}.haplotypecaller.vcf --reference $REF -ERC GVCF
gatk VariantFiltration --java-options "-Xmx8g" --variant $DIR/${PREFIX}.haplotypecaller.vcf --output $DIR/${PREFIX}.ASE.final.vcf --cluster-window-size 35 --cluster-size 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" 2>$DIR/$PREFIX.VariantFilter.log
