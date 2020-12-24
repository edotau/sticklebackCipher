#!/bin/sh
#SBATCH --mem=24G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --job-name=fqToBam
#SBATCH --exclude=dl-01
set -e

READ1=$1
READ2=$2

#Set reference:
REF=$3
export PATH=/data/lowelab/edotau/bin/envs/htslib/bin:$PATH
PREFIX=$(echo $READ1 | rev | cut -d '_' -f 2- | rev)
####Final output####

DIR=bamfiles
mkdir -p $DIR
BAM=$DIR/${PREFIX}.bam

#MAPPING

SORT_THREADS=$(($SLURM_CPUS_ON_NODE / 4))

THREADS=$(($SLURM_CPUS_ON_NODE - $SORT_THREADS))
SORTED_BAMS=$DIR/${PREFIX}.sorted.bam
echo -e "["$(date)"]\tAligning and sorting..."
echo "
bwa mem -t $THREADS $REF $READ1 $READ2 | samtools sort -@ $SORT_THREADS -T $DIR/${PREFIX}.tmp -m 1G /dev/stdin > $SORTED_BAMS"
bwa mem -t $THREADS $REF $READ1 $READ2 | samtools sort -@ $SORT_THREADS -T $DIR/${PREFIX}.tmp -m 1G /dev/stdin > $SORTED_BAMS

echo "
samtools index -@ $SLURM_CPUS_ON_NODE $SORTED_BAMS"
samtools index -@ $SLURM_CPUS_ON_NODE $SORTED_BAMS

module load GATK/4.1.3.0-gcb01
markedDups=$DIR/${PREFIX}.markedDups.bam

echo "
gatk MarkDuplicates -I $SORTED_BAMS -O $markedDups -M /dev/null -VALIDATION_STRINGENCY SILENT -AS"
gatk MarkDuplicates -I $SORTED_BAMS -O $markedDups -M /dev/null -VALIDATION_STRINGENCY SILENT -AS
samtools index -@ $SORT_THREADS $markedDups

library=$(echo $PREFIX | cut -d '_' -f 2 | cut -d '.' -f 1)
unit=$(samtools view $BAM | head -n 1 | cut -f 1 | cut -d ':' -f 3,4 | tr ':' '.')
platform="HiSeqX"
sampleID=$PREFIX

echo "
gatk AddOrReplaceReadGroups -I $markedDups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $sampleID | gatk FixMateInformation -I /dev/stdin -O $BAM -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR $DIR/"fixmate.tmp" -AS"
gatk AddOrReplaceReadGroups -I $markedDups -O /dev/stdout -LB $library -PL $platform -PU $unit -SM $sampleID | gatk FixMateInformation -I /dev/stdin -O $BAM -ADD_MATE_CIGAR true -IGNORE_MISSING_MATES true -TMP_DIR $DIR/"fixmate.tmp" -AS

samtools index -@ $SLURM_CPUS_ON_NODE $output
