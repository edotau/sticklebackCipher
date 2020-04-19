#!/bin/sh -e
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --job-name=toGenotypeVcf
#Set reference:
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta

#set number of threads to use, needs to equal computer or number you request in slurm,
#I found 8-12 usually works for me, anymore than 12 seems to slow things down
threads=12
module load bwa samtools/1.9-gcb01 GATK/4.1.3.0-gcb01
#used by trim_galore
module load cutadapt/2.3-gcb01 python/3.7.4-gcb01 pigz/2.3.4-gcb01

READ1=$1
READ2=$2

PREFIX=$(basename $READ1 _1.fq.gz )
DIR=${PREFIX}"_genotypeVcfs"
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
bwa mem -t $threads $REF $1 $2 | samtools view -hb /dev/stdin | samtools sort -@ $threads -T $DIR/${PREFIX}.tmp -o $BAM
samtools index $BAM

##Convert mapped bam into a bam that is a valid input for GATK
markedDups=$DIR/${PREFIX}.markedDups.bam
fixedBAM=$DIR/${PREFIX}.gatk.valid.bam
gVcf=$DIR/${PREFIX}.wgs.g.vcf.gz

gatk --java-options "-Xmx16G" MarkDuplicates -I $BAM -O $markedDups -M $DIR/${PREFIX}".dup.metrics" -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true
#adds Read group information to bam alignments
#this is how GATK differentiates between your cohort of samples
#required aruguments:
#Set real library read group if you need this field, for example, rna-seq lib, whole genome library atac, etc
library="Stickleback"
platform="Illumina"
unit="HiSeqX"
sampleID=$PREFIX

gatk --java-options "-Xmx16G" AddOrReplaceReadGroups -I $markedDups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $PREFIX | gatk --java-options "-Xmx16G" FixMateInformation -I /dev/stdin -O $fixedBAM -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true
samtools index $fixedBAM
#Now run GATK on processed BAM
gatk HaplotypeCaller --java-options "-Xmx32G" \
	-I $fixedBAM -O $gVcf -R $REF \
	-ERC GVCF --native-pair-hmm-threads $threads \
