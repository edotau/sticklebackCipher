#!/bin/sh
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --exclude=dl-01
set -e
target=$1
query=$2
PREFIX=$(basename $target .fa)_$(basename $query .fa)
axt=${PREFIX}.axt

echo "Target=$target, Query=$query"
export PATH=/data/lowelab/edotau/kentUtils/:$PATH

lastz=/data/lowelab/edotau/software/lastz-distrib-1.04.00/src/lastz
humanChimp=/data/lowelab/edotau/software/lastz-distrib-1.04.00/humanChimpMatrix.txt

echo "
$lastz $target $query --format=axt --output=$axt --scores=$humanChimp O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
"
$lastz $target $query --format=axt --output=$axt --scores=$humanChimp O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0 

echo "finished lastz"

echo "
axtChain -linearGap=medium -scoreScheme=$humanChimp $axt -faT $target -faQ $query /dev/stdout | chainSort /dev/stdin ${PREFIX}.chain
"
axtChain -linearGap=medium -scoreScheme=$humanChimp $axt -faT $target -faQ $query /dev/stdout | chainSort /dev/stdin ${PREFIX}.chain

echo "DONE"
