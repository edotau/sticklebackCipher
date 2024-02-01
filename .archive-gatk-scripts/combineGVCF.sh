#!/bin/sh
#SBATCH --mem=32G
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --job-name=MATAxRABS
#SBATCH --time=2-0

module load GATK/4.1.3.0-gcb01
PREFIX=$1
DIR=${2%/}
#Write to file, then remove it to make sure it did not contain samples from a different run.
samples=${PREFIX}.list
echo $PREFIX >> $samples
rm $samples
for i in $DIR/*.g.vcf.gz
do
	echo $i >> $samples
done
merged=${PREFIX}.merged.g.vcf.gz
VCF=${PREFIX}.cohort.GENOTYPED.vcf.gz
gatk CombineGVCFs --java-options "-Xmx30g" \
	-R $3 \
	-V $samples \
	-O $merged
rm $samples
gatk --java-options "-Xmx30g" GenotypeGVCFs \
	-R $3 \
	-V $merged \
	-O $VCF
#sbatch /data/lowelab/edotau/sticklebackCipher/gatk/getSNPS.sh $VCF
#sbatch /data/lowelab/edotau/sticklebackCipher/gatk/getIndels.sh $VCF
