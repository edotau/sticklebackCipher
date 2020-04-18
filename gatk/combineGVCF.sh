#!/bin/sh
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --job-name=mergeGenotypeVcf

module load GATK/4.1.3.0-gcb01
PREFIX=$1
#Write to file, then remove it to make sure it did not contain samples from a different run.
samples=${PREFIX}.list
echo $PREFIX >> $samples
rm $samples
for i in *.g.vcf.gz
do
	echo $i >> $samples
done
gatk CombineGVCFs --java-options "-Xmx8g" \
	-R /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta \
	-V $samples \
	-O ${PREFIX}.merged.vcf.gz
rm $samples