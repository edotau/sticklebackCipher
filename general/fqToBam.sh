#!/bin/sh
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --time=0-11:00:00
set -e
#Set reference:
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
module load bwa samtools/1.9-gcb01
#used by trim_galore
module load cutadapt/2.3-gcb01 python/3.7.4-gcb01 pigz/2.3.4-gcb01

READ1=$1
READ2=$2
#depending on which sequencing facility you use, fastq files might come in different naming conventions: ie read.fastq.gz or read.
#ex. echo CL12w16-3_atac_R1.fastq.gz | rev | cut -d '_' -f 2- | rev > CL12w16-3_atac
PREFIX=$(echo $READ1 | rev | cut -d '_' -f 2- | rev)
####Final output####
DIR="basic_alignment"
BAM=$DIR/${PREFIX}.bam
mkdir -p $DIR
fastqDIR="QCedFastqs"
trim_galore=/data/lowelab/edotau/software/TrimGalore-0.6.5/trim_galore
$trim_galore -o $fastqDIR --basename $PREFIX --trim-n --max_n 0 --cores 4 --paired $READ1 $READ2
trim1=$fastqDIR/${PREFIX}_val_1.fq.gz
trim2=$fastqDIR/${PREFIX}_val_2.fq.gz

input1=$fastqDIR/${PREFIX}_R1.fastq.gz
input2=$fastqDIR/${PREFIX}_R2.fastq.gz

#Second Trimming of adaptors:
/data/lowelab/software/bbmap/bbduk.sh \
	in1=$trim1 in2=$trim2 \
	out1=$input1 out2=$input2 \
	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
	ref=/data/lowelab/software/bbmap/resources/adapters.fa

#MAPPING
output=$DIR/${PREFIX}.FINAL.bam

SORT_THREADS=$(($SLURM_CPUS_ON_NODE / 4))

THREADS=$(($SLURM_CPUS_ON_NODE - $SORT_THREADS))


echo -e "["$(date)"]\tAligning and sorting..."
echo "
bwa mem -t $THREADS $REF $READ1 $READ2 | samtools sort -@ $SORT_THREADS -T $DIR/${PREFIX}.tmp -m 1G /dev/stdin > $BAM"
bwa mem -t $THREADS $REF $READ1 $READ2 | samtools sort -@ $SORT_THREADS -T $DIR/${PREFIX}.tmp -m 1G /dev/stdin > $BAM

echo "
samtools index -@ $SLURM_CPUS_ON_NODE $BAM"
samtools index -@ $SLURM_CPUS_ON_NODE $BAM

module load GATK/4.1.3.0-gcb01
markedDups=$DIR/${PREFIX}.markedDups.bam

echo "
gatk MarkDuplicates -I $BAM -O $markedDups -M /dev/null -VALIDATION_STRINGENCY SILENT -AS"
gatk MarkDuplicates -I $BAM -O $markedDups -M /dev/null -VALIDATION_STRINGENCY SILENT -AS
samtools index -@ $SORT_THREADS $markedDups

library=$(echo $PREFIX | cut -d '_' -f 2 | cut -d '.' -f 1)
unit=$(samtools view $BAM | head -n 1 | cut -f 1 | cut -d ':' -f 3,4 | tr ':' '.')
platform="HiSeqX"
sampleID=$PREFIX


echo "
gatk AddOrReplaceReadGroups -I $markedDups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $sampleID | gatk FixMateInformation -I /dev/stdin -O $output -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR $DIR/"fixmate.tmp" -AS"
gatk AddOrReplaceReadGroups -I $markedDups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $sampleID | gatk FixMateInformation -I /dev/stdin -O $output -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR $DIR/"fixmate.tmp" -AS

samtools index -@ $SLURM_CPUS_ON_NODE $output

exit 0
