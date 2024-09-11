#!/bin/bash -e

# Run lastz
DIR="testdata"

bash lastz.sh $DIR/freshwater.fa $DIR/marine.fa

bash mergeSortChains.sh $DIR/freshwater.fa $DIR/marine.fa ./

bash liftover.sh $DIR/target.fa $DIR/query.fa freshwater.marine.all.sorted.chain

bash annotations.sh $DIR/freshwater.fa freshwater.marine.all.sorted.chain $DIR/chrM.gtf  output