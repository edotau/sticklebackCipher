#!/bin/sh -e
module add GATK/4.1.3.0-gcb01
gatk --java-options "-Xmx16g" GenotypeGVCFs -R /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.fasta -V $1 -O $(basename $1 .g.vcf)GenotypeGVCFs.vcf
