#!/bin/sh
#SBATCH --mem=16G 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eric.au@duke.edu
#SBATCH --parsable
#SBATCH --job-name=hic_bwa_mapping01
module load samtools bwa perl/5.10.1-fasrc04 java/1.8.0_45-fasrc01
FILTER='/data/lowelab/edotau/software/mapping_pipeline/filter_five_end.pl'
input=$1
REF=$2
PREFIX=$(echo $input | sed 's/.fastq.gz//')
bwa mem -t 12 -B 8 $REF $input | samtools view -Sb - > ${PREFIX}.bam
samtools view -h ${PREFIX}.bam | perl $FILTER | samtools view -Sb - > ${PREFIX}_filter.bam
