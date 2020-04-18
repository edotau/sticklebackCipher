#!/bin/sh
module load GATK/4.1.0.0-gcb01 java/1.8.0_45-fasrc01 samtools/1.9-gcb01 fastqc/0.11.5-fasrc01 STAR/2.7.2b-gcb01
picard=/data/lowelab/edotau/software/picard.jar

#RGPU=EA02
#library=atac

#PREFIX=$(echo $i | sed 's/.bam//')

#star_ref=/data/lowelab/RefGenomes/gasAcu1/STARindex
#ref=/data/lowelab/RefGenomes/gasAcu1/gasAcu1.fa
star_ref=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/STARidx

REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
input1=$1
input2=$2

PREFIX=$(echo $input1 | sed 's/w1_L8_R1_val_1.fq.gz//')toRABS
bn=$PREFIX
opdir=$bn"rnaseq"
mkdir $opdir

READ1=$opdir/${PREFIX}_R1.fastq.gz
READ2=$opdir/${PREFIX}_R2.fastq.gz

#fastqc -t 6 $READ1 $READ2 -o $opdir
#Preprocessing
#Trimming adaptors
/data/lowelab/software/bbmap/bbduk.sh \
	in1=$input1 \
	in2=$input2 \
	out1=$READ1 \
	out2=$READ2 \
	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
	ref=/data/lowelab/software/bbmap/resources/adapters.fa

#MAPPING

echo -e "["$(date)"]\tAligning.."
STAR --outFileNamePrefix $opdir/$bn --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMattrRGline ID:$bn CN:Gonomics LB:PairedEnd PL:Illumina PU:Unknown SM:$bn --genomeDir $star_ref --runThreadN 8 --readFilesCommand zcat --readFilesIn $READ1 $READ2 --twopassMode Basic

echo -e "["$(date)"]\tSorting.."
samtools sort -o $opdir/$bn"_sorted.bam" $opdir/$bn"Aligned.out.bam"
rm $opdir/$bn"Aligned.out.bam"

echo -e "["$(date)"]\tIndexing.."
samtools index $opdir/$bn"_sorted.bam"

echo -e "["$(date)"]\tMarking duplicates.."
java -Xmx8G -jar $picard MarkDuplicates I=$opdir/$bn"_sorted.bam" O=$opdir/$bn"_dupMarked.bam" M=$opdir/$bn"_dup.metrics" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT 2>$opdir/$bn.MarkDuplicates.log
rm $opdir/$bn"_sorted.bam"
rm $opdir/$bn"_sorted.bam.bai"

#SplitNCigarReads
echo -e "["$(date)"]\tSpliting reads.."
gatk SplitNCigarReads --java-options "-Xmx8g" --input $opdir/$bn"_dupMarked.bam" --output $opdir/$bn"_split.bam" --reference $REF 2>$opdir/$bn.SplitNCigarReads.log
samtools index $opdir/$bn"_split.bam"


gatk HaplotypeCaller --java-options "-Xmx8g" --input $opdir/$bn"_split.bam" --dont-use-soft-clipped-bases true --output $opdir/${bn}.haplotypecaller.vcf --reference $REF -ERC GVCF
gatk VariantFiltration --java-options "-Xmx8g" --variant $opdir/${bn}.haplotypecaller.vcf --output $opdir/${bn}.ASE.final.vcf --cluster-window-size 35 --cluster-size 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" 2>$opdir/$bn.VariantFilter.log
