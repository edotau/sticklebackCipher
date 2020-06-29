#!/bin/bash
REF=$1
BAM=$2
sbatch ./gatkHaplotype.sh $BAM $REF

#!/bin/bash
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4 --ntasks=3
#SBATCH --array=1-100
#SBATCH --time=0-03
#SBATCH --nodes=1
#SBATCH --exclude=dl-01,x2-02-2,c1-10-3

i=$SLURM_ARRAY_TASK_ID
bam=$1
fa=$2
faidx=$fa.fai
#LEN=$3

#####tasks=$SLURM_CPUS_ON_NODE
LEN="$((wc -l $contigs))"

# Perform function every %100 = $i th line
for j in $(seq $i 100 $LEN )
do
    echo $j
    contig=`sed -n ${j}p $faidx | awk '{print $1}'`
    contig_no_pipe=`echo $contig | sed 's/|/_/g'`
    end=`sed -n ${j}p $faidx | awk '{print $2}'`
    
    #####LEN="$((wc -l $faidx)+$tasks-1/$tasks))"

    chunk=$(($LEN / 3))
    #last chuck includes remainder
    s2=$start+$chunk, s3=$s2+$chunk , s4=$s3+$end
    gatk $bam -L $contig:$start-$s2 & gatk $bam -L $contig:$s2-$s3 & gatk $bam -L $contig:$s3-$end
done
exit 0
