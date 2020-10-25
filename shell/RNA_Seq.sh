#!/bin/sh
#SBATCH --mem=16G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --time=2-0
#SBATCH --job-name=Snps.RNA-Seq
#SBATCH --exclude=dl-01
set -e
module load GATK/4.1.0.0-gcb01 java/1.8.0_45-fasrc01 samtools/1.9-gcb01 fastqc/0.11.5-fasrc01 STAR/2.7.2b-gcb01
module load cutadapt/2.3-gcb01 python/3.7.4-gcb01 pigz/2.3.4-gcb01
picard=/data/lowelab/edotau/software/picard.jar

#RGPU=EA02
#library=atac

PREFIX=$(echo $i | sed 's/.bam//')

star_ref=/data/lowelab/edotau/toGasAcu1.5/idx/STAR_2

REF=/data/lowelab/edotau/toGasAcu1.5/idx/stickleback_v5_assembly.fa
input1=$1
input2=$2

PREFIX=$(echo $input1 | cut -d '_' -f 1)
bn=$PREFIX
opdir=$bn".RNA-Seq"
mkdir -p $opdir

fastqc -t $SLURM_CPUS_ON_NODE $input1 $input2 -o $opdir
#Preprocessing
#Trimming adaptors using kmer hash
firstTrimOne=${PREFIX}_R1.fastq.gz
firstTrimTwo=${PREFIX}_R2.fastq.gz
/data/lowelab/edotau/bin/github.com/bbmap/bbduk.sh \
	in1=$input1 \
	in2=$input2 \
	out1=$firstTrimOne \
	out2=$firstTrimTwo \
	minlen=25 qtrim=rl trimq=10 ktrim=r k=25 mink=11 hdist=1 \
	ref=/data/lowelab/edotau/bin/github.com/bbmap/resources/adapters.fa
	
fastqDIR="QCedFastqs"
#Trimming adaptors remaining from kmer hash, wraps cutadapt
$trim_galore -o $fastqDIR --basename $PREFIX --trim-n --max_n 0 --cores 4 --paired $firstTrimOne $firstTrimTwo
trim1=$fastqDIR/${PREFIX}_val_1.fq.gz
trim2=$fastqDIR/${PREFIX}_val_2.fq.gz
#MAPPING splice aware aligner

echo -e "["$(date)"]\tAligning.."
STAR --outFileNamePrefix $opdir/$bn --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMattrRGline ID:$bn CN:Gonomics LB:PairedEnd PL:Illumina PU:HiSeqX SM:$bn --genomeDir $star_ref --runThreadN $SLURM_CPUS_ON_NODE --readFilesCommand zcat --readFilesIn trim1 trim2 --twopassMode Basic

echo -e "["$(date)"]\tSorting.."
samtools sort -@ $SLURM_CPUS_ON_NODE -o $opdir/$bn"_sorted.bam" $opdir/$bn"Aligned.out.bam"
rm $opdir/$bn"Aligned.out.bam"

echo -e "["$(date)"]\tIndexing.."
samtools index -@ $SLURM_CPUS_ON_NODE $opdir/$bn"_sorted.bam"

echo -e "["$(date)"]\tMarking duplicates.."
java -jar $picard MarkDuplicates I=$opdir/$bn"_sorted.bam" O=$opdir/$bn"_dupMarked.bam" M=$opdir/$bn"_dup.metrics" VALIDATION_STRINGENCY=SILENT 2>$opdir/$bn.MarkDuplicates.log
rm $opdir/$bn"_sorted.bam"
rm $opdir/$bn"_sorted.bam.bai"

samtools index -@ $SLURM_CPUS_ON_NODE $opdir/$bn"_dupMarked.bam"

#SplitNCigarReads
echo -e "["$(date)"]\tSpliting reads.."
gatk SplitNCigarReads --input $opdir/$bn"_dupMarked.bam" --output $opdir/$bn"_split.bam" --reference $REF 2>$opdir/$bn.SplitNCigarReads.log
samtools index -@ $SLURM_CPUS_ON_NODE $opdir/$bn"_split.bam"

gatk HaplotypeCaller --input $opdir/${bn}_split.bam --dont-use-soft-clipped-bases true --output $output --reference $REF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation -ERC GVCF

echo "RNA transcript ASSEMBLY...(stringtie)"

module add StringTie/2.1.1-gcb01
gtf=/data/lowelab/edotau/rabsTHREEspine/index/bepa_v1_ensembl_pslMappedRabs.gtf
#/data/lowelab/edotau/toGasAcu1.5/idx/stickleback_v5_maker_genes_n
stringtie $opdir/$bn"_dupMarked.bam" -p $SLURM_CPUS_ON_NODE -G $gtf -o ${PREFIX}.RNA-Seq.gtf -l ${PREFIX^^}


#submit=${PREFIX}.gatk.txt
#/data/lowelab/edotau/toGasAcu1.5/scripts/gatkToText.sh $opdir/$bn"_split.bam" $REF $submit
#/data/lowelab/edotau/toGasAcu1.5/scripts/array.gatk.sh $submit


exit 0


#gatk VariantFiltration --java-options "-Xmx8g" --variant $opdir/${PREFIX}.haplotypecaller.vcf --output $opdir/${bn}.ASE.final.vcf --cluster-window-size 35 --cluster-size 3 --filter-name FS --filter-expression "FS > 30.0" --filter-name QD --filter-expression "QD < 2.0" 2>$opdir/$bn.VariantFilter.log
