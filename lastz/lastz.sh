#!/bin/sh
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
set -e
# Set Path to kentUtils or load kentUtils module if available on your cluster:
export PATH=/data/lowelab/edotau/kentUtils/:$PATH

# Set lastz or add to path
lastz=/data/lowelab/edotau/software/lastz-distrib-1.04.00/src/lastz

# Determine Scoring matrix
scoreMatrix=/data/lowelab/edotau/software/lastz-distrib-1.04.00/humanChimpMatrix.txt

# target or reference fasta MUST be a single record fasta sequence:
target=$1
# query fasta can be multi records
query=$2  # fasta query

if [[ "$#" -lt 2 ]]; then
    echo "Usage: ./lastz.sh target.fa query.fa"
    exit 0
fi

PREFIX=$(basename $target .fa)_$(basename $query .fa)
axt=${PREFIX}.axt

echo "
$lastz $target $query --format=axt --output=$axt --scores=$scoreMatrix O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0
"
$lastz $target $query --format=axt --output=$axt --scores=$scoreMatrix O=600 E=150 T=2 M=254 K=4500 L=4500 Y=15000 C=0 

echo "finished lastz"

echo "
axtChain -linearGap=medium -scoreScheme=$scoreMatrix $axt -faT $target -faQ $query /dev/stdout | chainSort /dev/stdin ${PREFIX}.chain
"
axtChain -linearGap=medium -scoreScheme=$scoreMatrix $axt -faT $target -faQ $query /dev/stdout | chainSort /dev/stdin ${PREFIX}.chain

echo "DONE"
