#!/bin/sh
canu=/data/lowelab/edotau/software/canu/Linux-amd64/bin/canu
module add java/9-gcb01
module add gnuplot/4.6.4-fasrc02
script=/data/lowelab/edotau/data/pacbio_rabs_688/run_canu.sh
sbatch \
	--mem=100G \
	--nodes=1 \
	--cpus-per-task=8 \
	--ntasks=2 \
	--mem-per-cpu=8g \
	--wrap="$canu -p RS688.canu \
		-d RS688_nanoPacBio.canu \
		corThreads=8 \
		corConcurrency=16 \
		minReadLength=10000 \
		minOverlapLength=1000 \
		corMhapSensitivity=high \
		utgmhapThreads=16 \
		batThreads=16 \
		cnsThreads=16 \
		useGrid=1 \
		gridEngine=slurm \
		gridOptionsUTGMHAP=16 \
		gridOptionsutgovl=16 \
		gridOptionsCNS=16 \
		gridOptionsBAT=16 \
		gridOptionsUTGMHAP=16 \
		gridOptionsExecutive=--mem=64g \
		genomeSize=600m \
		overlapper=mhap \
		utgReAlign=true \
		gridOptions="--cpus-per-task=8 --ntasks=3 --mem-per-cpu=8g" \
		mhapThreads=16 \
		-pacbio-raw /data/lowelab/edotau/data/pacbio_rabs_688/raw_bam_alignment/*.fastq* \
		-nanopore-raw /data/lowelab/edotau/stickleback_genome_raw_data/20181211_0015_RABS688_Nanopore/*.fastq*"
