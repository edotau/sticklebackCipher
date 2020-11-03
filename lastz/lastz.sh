#!/bin/sh
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
target=$1
query=$2
module add kentUtils/v302-gcb01
lastz=/data/lowelab/edotau/software/lastz-distrib-1.04.00/src/lastz
humanChimp=/data/lowelab/edotau/software/lastz-distrib-1.04.00/humanChimpMatrix.txt
PREFIX=$(basename $target .fa)_$(echo $query | cut -d '.' -f 1)
echo $target
axt=${PREFIX}.axt
$lastz $target $query --format=axt --output=$axt --scores=$humanChimp \
        O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
echo "finished lastz"
axtChain -linearGap=medium -scoreScheme=$humanChimp $axt -faT $target -faQ $query /dev/stdout | chainSort /dev/stdin $chainOut ${PREFIX}.chain
echo "DONE"
