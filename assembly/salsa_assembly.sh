#!/bin/sh
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --parsable
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eric.au@duke.edu
#SBATCH --job-name=HicAssembly
module load python/2.7.11-fasrc01 boost/1.59.0-gcb01 gcc/6.2.0-fasrc01 java samtools
salsa=/data/lowelab/edotau/sticklebackCipher/assembly/SALSA/run_pipeline.py

bedfile=$1
REF=$2
FAIDX=$REF.fai
samtools faidx $REF
python $salsa -a $REF -l $FAIDX -b $bedfile -e GATC,GANTC -m yes -p yes -c 10000 -o $3

echo "Finished scaffolding, now generating hi-c contact maps, could take several days...."
/data/lowelab/edotau/software/SALSA/convert.sh $3/

#stitch=/data/lowelab/edotau/software/SALSA/stitch.py
#unitig_bed=/data/lowelab/edotau/scratch/arimaHiC/unitigs_canu/rabs.4.1.canu.unitigs.bed
#unitig_file=/data/lowelab/edotau/scratch/arimaHiC/unitigs_canu/rabs.4.1.canu.unitigs.fasta

#sbatch --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=8 --mail-type=ALL --mail-user=eric.au@duke.edu /data/lowelab/edotau/scratch/arimaHiC/unitigs_canu/fastq/stitch.sh
#$stitch -c unitig_pacbio_nano/scaffolds_FINAL.fasta -b $unitig_bed -u $unitig_file -p yes -o unitig_pacbio_nanoFINAL
#/data/lowelab/edotau/software/SALSA/convert.sh unitig_pacbio_nano

#-u /data/lowelab/edotau/scratch/arimaHiC/canu_pabBioNano/ref/rabs.4.1.canu.unitigs.bed \
#	-g /data/lowelab/edotau/scratch/arimaHiC/unitigs_canu/rabs.4.1.canu.unitigs.gfa
