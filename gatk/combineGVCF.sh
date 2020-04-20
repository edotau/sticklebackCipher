#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --job-name=mergeGvcfToGenotype

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
merged=${PREFIX}.merged.g.vcf.gz
VCF=${PREFIX}.GENOTYPED.merged.vcf.gz
gatk CombineGVCFs --java-options "-Xmx32g" \
	-R /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta \
	-V $samples \
	-O $merged
rm $samples
gatk --java-options "-Xmx30g" GenotypeGVCFs \
	-R /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta \
	-V $merged \
	-O $VCF
sbatch /data/lowelab/edotau/sticklebackCipher/gatk/getSNPS.sh $VCF
sbatch /data/lowelab/edotau/sticklebackCipher/gatk/getIndels.sh $VCF
