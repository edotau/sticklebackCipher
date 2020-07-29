#!/bin/sh
#SBATCH --mem=50G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --time=2-12
#SBATCH --job-name=litcMata
#SBATCH --exclude=dl-01,c1-10-3
set -e
#Set reference:

READ1=$1
READ2=$2
REF=$3
#depending on which sequencing facility you use, fastq files might come in different naming conventions: ie read.fastq.gz or read.
#ex. echo CL12w16-3_atac_R1.fastq.gz | rev | cut -d '_' -f 2- | rev > CL12w16-3_atac
PREFIX=$(echo $READ1 | rev | cut -d '_' -f 2- | rev)
####Final output####
#bash /data/lowelab/edotau/bin/etc/profile.d/conda.sh
DIR=$4
mkdir -p $DIR
BAM=$DIR/${PREFIX}.bam
#fastqDIR="QCedFastqs"
#trim_galore=/data/lowelab/edotau/software/TrimGalore-0.6.5/trim_galore
#$trim_galore -o $fastqDIR --basename $PREFIX --trim-n --max_n 0 --cores 4 --paired $READ1 $READ2
#trim1=$fastqDIR/${PREFIX}_val_1.fq.gz
#trim2=$fastqDIR/${PREFIX}_val_2.fq.gz

#input1=$fastqDIR/${PREFIX}_R1.fastq.gz
#input2=$fastqDIR/${PREFIX}_R2.fastq.gz

#Second Trimming of adaptors:
#/data/lowelab/software/bbmap/bbduk.sh \
#	in1=$trim1 in2=$trim2 \
#	out1=$input1 out2=$input2 \
#	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
#	ref=/data/lowelab/software/bbmap/resources/adapters.fa

#clean directory: rm $trim1 $ trim2
#MAPPING
echo -e "["$(date)"]\tAligning and sorting..."
bwa mem -t $SLURM_CPUS_ON_NODE $REF $READ1 $READ2 | samtools sort -@ $SLURM_CPUS_ON_NODE -T $DIR/${PREFIX}.tmp -m 2G /dev/stdin > $BAM
samtools index -@ $SLURM_CPUS_ON_NODE $BAM
module load GATK/4.1.3.0-gcb01
markedDups=$DIR/${PREFIX}.markedDups.bam
output=$5
#$DIR/${PREFIX}.gatk.valid.bam
gatk --java-options "-Xmx32G" MarkDuplicates -I $BAM -O $markedDups -M /dev/null -VALIDATION_STRINGENCY SILENT -CREATE_INDEX true

library=$(echo $PREFIX | cut -d '_' -f 2 | cut -d '.' -f 1)
unit=$(samtools view $BAM | head -n 1 | cut -f 1 | cut -d ':' -f 3,4 | tr ':' '.')
platform="HiSeqX"
sampleID=$PREFIX

gatk --java-options "-Xmx32G" AddOrReplaceReadGroups -I $markedDups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $sampleID | gatk --java-options "-Xmx32G" FixMateInformation -I /dev/stdin -O $output -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR $DIR/"fixmate.tmp"

samtools index -@ $SLURM_CPUS_ON_NODE $output

/data/lowelab/edotau/RABSGenome/gatk/haplotype.sh $output $REF

exit 0
