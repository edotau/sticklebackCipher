#!/bin/sh -e
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --job-name=toGenotypeVcf
#SBATCH --parsable
BAM=$1
PREFIX=$(basename $BAM .bam)

#Either load modules or assign paths to GATK, samtools, and picard tools
module load GATK/4.1.3.0-gcb01
#Set your reference:
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
#if you want to name your files anything else set these variables

DIR=$2
"toGenotypeVcfs"
mkdir -p $DIR
markedDups=$DIR/${PREFIX}.markedDups.bam
output=$DIR/${PREFIX}.gatk.valid.bam
#final output
gVcf=${PREFIX}.g.vcf.gz

#Now run GATK on processed BAM
gatk HaplotypeCaller --java-options "-Xmx32G" \
	-I $output -O $gVcf -R $REF \
	-ERC GVCF --native-pair-hmm-threads 12 \
rm $markedDups ${markedDups}.bai
