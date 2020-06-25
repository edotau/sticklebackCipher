#!/bin/bash
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4 --ntasks=2
#SBATCH --array=1-100
#SBATCH --time=0-03
#SBATCH --nodes=1
#SBATCH --exclude=dl-01,x2-02-2,x2-04-4

i=$SLURM_ARRAY_TASK_ID

bam=$1

LEN=`wc -l $fasta.fai | awk '{print $1}'`

# Perform function every %100 = $i th line
for j in $(seq $i 100 $LEN )
do
    echo $j
    contig=`sed -n ${j}p $fasta.fai | awk '{print $1}'`
    contig_no_pipe=`echo $contig | sed 's/|/_/g'`
    end=`sed -n ${j}p $fasta.fai | awk '{print $2}'`
    
    #chunk=$end%4
    s2=$start+$chunk, s3=$s2+$chunk , s4=$s3+$chunk
    gatk $bam -L $contig:$start-$s2 & gatk $bam -L $contig:$s2-$s3 & gatk $bam -L $contig:$s3-$end
done
