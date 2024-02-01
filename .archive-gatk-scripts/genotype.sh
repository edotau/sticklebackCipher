#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --job-name=getGenotypes
module load GATK/4.1.3.0-gcb01
gVCF=$1
PREFIX=$(basename $gVCF .g.vcf.gz)
REF=/data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta
gatk --java-options "-Xmx30g" GenotypeGVCFs -R $REF -V $gVCF -O ${PREFIX}.vcf.gz
