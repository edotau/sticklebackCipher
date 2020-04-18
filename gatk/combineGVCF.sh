#!/bin/sh
module load GATK/4.1.3.0-gcb01
PREFIX=$1
#Write to file, then remove it to make sure it did not contain samples from a different run.
samples=${PREFIX}.list
echo $PREFIX >> $samples
rm $samples
for i in *.g.vcf
do
	echo $i >> $samples
done
gatk CombineGVCFs --java-options "-Xmx8g" \
	-R /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta \
	-V $samples \
	-O ${PREFIX}.merged.vcf.gz
