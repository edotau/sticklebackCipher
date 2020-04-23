#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --job-name=getGenotypes
module load GATK/4.1.3.0-gcb01
BAM=$1
PREFIX=$(basename $BAM .bam)
gatk HaplotypeCaller --java-options "-Xmx32G" --input $BAM --output ${PREFIX}.g.vcf --reference /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta -ERC GVCF
