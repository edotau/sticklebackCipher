#!/bin/sh -e
bed=$1
module add gcc kentUtils/v302-gcb01
PREFIX=$(basename $bed .bed)LiftoverToRabs
liftoverChain=/data/lowelab/edotau/toGasAcu2RABS/annotation/gasAcu1ToRabsLiftover.chain
#-minMatch=0.85
liftOver $bed $liftoverChain ${PREFIX}Temp.bed ${PREFIX}UnmappedToRabs
sort -k1,1 -k2,2n ${PREFIX}Temp.bed > ${PREFIX}.bed
bedToBigBed ${PREFIX}.bed /data/lowelab/edotau/toGasAcu2RABS/gasAcu2RABS/gasAcu2RABS.sizes ${PREFIX}.bb
rm ${PREFIX}Temp.bed
scp ${PREFIX}.bb eha17@trackhub.genome.duke.edu:/nfs/trackhub/lowelab/edotau/myHub/annotation
