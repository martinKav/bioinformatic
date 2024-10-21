#!/bin/bash
projPath="/mnt/c/users/mjmka/Documents/Bioinformatics" 
cores=8
ref="/mnt/c/users/mjmka/Documents/Bioinformatics/bowtie2Index/hg38/GRCh38_noalt_as"
picardCMD="java -jar $projPath/tools/picard/build/libs/picard.jar"
binLen=500
chromSize="/mnt/c/users/mjmka/Documents/Bioinformatics/ncbi_dataset/hg38_genome_size.txt"
seacr="/mnt/c/users/mjmka/Documents/Bioinformatics/tools/SEACR/SEACR_1.3.sh" 
histList=("IPSC_H3K27ac" "IPSC_H3K27me3" "IPSC_yH2AX" "IPSC_2O" "D7_H3K27ac" "D7_H3K27me3" "D7_yH2AX" "D7_2O" "Undetermined")
 for histName in ${histList[@]} ; do
# rm -rf $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
# rm -f -rf $projPath/alignment/rmDuplicate/picard_summary/${histName}_picard.rmDup.txt
# rm -f -rf $projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.rmDup.txt
# rm -f -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragmentsCount.bin$binLen.bed
# rm -rf $projPath/alignment/sam/${histName}_bowtie2.sam
# rm -f -rf $projPath/alignment/sam/${histName}_bowtie2.sorted.sam
# rm -f -rf $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam
# rm -f -rf $projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam
# rm -f -rf $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam
# rm -f -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.bed
# rm -f -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed
# rm -f -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed

# rm -f -rf $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph


# rm -rf fastqFileQC/${histName}
# rm -rf fastqFileQC/${histName}_postcut
# rm -rf $projPath/fastq/${histName}_R1.QC.fastq.gz
# rm -rf $projPath/fastq/${histName}_R2.QC.fastq.gz
mkdir $projPath/alignments
mkdir $projPath/alignment/bam
mkdir $projPath/alignment/sam
mkdir $projPath/alignment/bed
mkdir $projPath/alignment/bedgraph
mkdir $projPath/fastq
mkdir $projPath/fastqFileQC
mkdir ${projPath}/alignment/sam/bowtie2_summary
mkdir $projPath/alignment/sam/fragmentLen
mkdir $projPath/alignment/rmDuplicate/picard_summary/
mkdir $projPath/alignment/rmDuplicate
mkdir $projPath/MK3/${histName}/
mkdir $projPath/fastqFileQC/${histName}
mkdir $projPath/fastqFileQC/${histName}_postcut

cp $projPath/data/${histName}/*_R1_001.gz $projPath/MK3/${histName}/${histName}_R1.fastq.gz
cp $projPath/data/${histName}/*_R2_001.gz $projPath/MK3/${histName}/${histName}_R2.fastq.gz



echo "cutadapting" ${histName}
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -j 8 -m 12 -q 25,25 --nextseq-trim 20 -o $projPath/fastq/${histName}_R1.QC.fastq.gz -p $projPath/fastq/${histName}_R2.QC.fastq.gz $projPath/MK3/${histName}/${histName}_R1.fastq.gz $projPath/MK3/${histName}/${histName}_R2.fastq.gz >${projPath}/fastq/cut_results/${histName}_cutadapt_summary.txt

echo "fastq examining cut" ${histName}
$projPath/tools/FastQC/fastqc -o ${projPath}/fastqFileQC/${histName}_postcut -f fastq ${projPath}/fastq/${histName}_R1.QC.fastq.gz
$projPath/tools/FastQC/fastqc -o ${projPath}/fastqFileQC/${histName}_postcut -f fastq ${projPath}/fastq/${histName}_R2.QC.fastq.gz

# not using bwa
# $projPath/tools/bwa/bwa mem -t 8 -M \
# -R "@RG\tID:${histName}\tPL:ILLUMINA\tSM:${histName}"\
# -P $projPath/bwaIndex/GCF_000001405.40_GRCh38.p14_genomic.fna  \
# $projPath/fastq/${histName}_R1.QC.fastq.gz $projPath/fastq/${histName}_R2.QC.fastq.gz \
# 2> $projPath/alignment/sam/bwa_log/${histName}_bwa.err \
# > $projPath/alignment/sam/${histName}_bwa.sam

echo "aligning with bowtie2" $histName
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-p ${cores} --rg-id ${histName} --rg SM:${histName} --rg PL:ILLUMINA --rg LB:MK3 --rg PU:MK3_${histName} -x ${ref} -1 \
${projPath}/fastq/${histName}_R1.QC.fastq.gz -2 ${projPath}/fastq/${histName}_R2.QC.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt

echo "sorting sam" $histName
$picardCMD SortSam -I $projPath/alignment/sam/${histName}_bowtie2.sam -O $projPath/alignment/sam/${histName}_bowtie2.sorted.sam --SORT_ORDER coordinate
echo "removing dupes" ${histName}
$picardCMD MarkDuplicates -I $projPath/alignment/sam/${histName}_bowtie2.sorted.sam -O $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam \
	--REMOVE_SEQUENCING_DUPLICATES true -M $projPath/alignment/rmDuplicate/picard_summary/${histName}_picard.rmDup.txt

echo "grabbing mapped size of dup-removed" $histName
samtools view -F 0x04 $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort |uniq -c |awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.rmDup.txt
echo "converting to bam" $histName
samtools view -h -bS -F 0x04 $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam >$projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam

 echo "removing intermiediate sams" ${histName}
 rm -rf $projPath/alignment/sam/${histName}_bowtie2.sam
 rm -rf $projPath/alignment/sam/${histName}_bowtie2.sorted.sam
 rm -rf $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam

echo "sorting bams" ${histName}
samtools sort -n $projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam -o $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam
echo "bam to bed" ${histName}
bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam -bedpe >$projPath/alignment/bed/${histName}_bowtie2.rmDup.bed
echo "cleaning bed" ${histName}
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.rmDup.bed >$projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed
echo "extracting frag columns" ${histName}
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed
echo "removing any leftover alt contigs" ${histName}
grep -E "chr1$(printf '\t')|chr2$(printf '\t')|chr3$(printf '\t')|chr4$(printf '\t')|chr5$(printf '\t')|chr6$(printf '\t')\
|chr7$(printf '\t')|chr8$(printf '\t')|chr9$(printf '\t')|chr10$(printf '\t')|chr11$(printf '\t')|chr12$(printf '\t')|chr13$(printf '\t')\
|chr14$(printf '\t')|chr15$(printf '\t')|chr16$(printf '\t')|chr17$(printf '\t')|chr18$(printf '\t')|chr19$(printf '\t')|chr20$(printf '\t')\
|chr21$(printf '\t')|chr22$(printf '\t')|chrX$(printf '\t')|chrY$(printf '\t')|chrM$(printf '\t')" $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed >$projPath/alignment/bed/${histName}_bowtie2.refiltered.bed
echo "assessing reproductibility"
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.refiltered.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.rmDup.fragmentsCount.bin$binLen.bed

echo "removing intermiediate bams and beds" ${histName}
rm $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam
rm -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.bed
rm -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed
rm -rf $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed
echo "done"
done

# echo "executing ChIPseqSpikeInFree script"
# rm -rf $projPath/alignment/ChIPseqSpikeInFree 
# mkdir $projPath/alignment/ChIPseqSpikeInFree
# Rscript $projPath/scripts/ChIPSpikeInFree.R $projPath
# rm -rf $projPath/alignment/ChIPseqSpikeInFree/SpikeInFree_results_rawCounts.txt

# for histName in ${histList[@]} ; do
# echo "getting scale factor" ${histName}
# scale_factor=`cat $projPath/alignment/ChIPseqSpikeInFree/${histName}.scaleFactor.txt`
# echo "bed to bedgraph" ${histName}
# bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.refiltered.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph
# echo "sorting bedgraph" ${histName}
# sort -k1,1 -k2,2n $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph >$projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.sorted.bedgraph
# echo "filtering bedgraph" ${histName}
# Rscript $projPath/scripts/removeBlacklist.R $projPath $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.sorted.bedgraph $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.badEOL.bedgraph 
# cat $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.badEOL.bedgraph | tr -d " \r" >$projPath/alignment/bedgraph/${histName}_bowtie2.filtered.bedgraph
# echo "sorting and merging bedgraph" ${histName}
# sort -k1,1V -k2,2n $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.bedgraph >$projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph
# bedtools merge -i $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph -c 4 -d 0 -o max > $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.bedgraph
# LC_COLLATE=C sort -k1,1 -k2,2n $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.bedgraph >$projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.sorted.bedgraph
# echo "converting bedgraph to bigwig" $histName
# bedGraphToBigWig $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.sorted.bedgraph $chromSize $projPath/alignment/bigWig/${histName}.rmDup.bigWig

# echo "removing intermiediate beds etc" ${histName}

# rm -rf $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph
# rm -rf $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.sorted.bedgraph
# rm -rf $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.badEOL.bedgraph 
# rm -rf $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.bedgraph
# rm -rf $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.bedgraph
# rm -rf $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.sorted.bedgraph
# done 

# for histControl in IPSC_2O D7_2O; do
# if [[ histControl == IPSC_2O ]]; then
# for histName in IPSC_H3K27ac IPSC_H3K27me3 IPSC_yH2AX D7_2O Undetermined; do
# echo "peakCalling for" histName "with" $histControl
# for peakNature in stringent relaxed; do
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph $projPath/alignment/bedgraph/${histControl}_bowtie2.filtered.sorted.bedgraph norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_control_${peakNature}.peaks
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_top0.01_${peakNature}.peaks
# done; done
# else
# for histName in IPSC_2O D7_H3K27ac D7_H3K27me3 D7_yH2AX;do
# echo "peakCalling for" histName "with" $histControl
# for peakNature in stringent relaxed; do
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph $projPath/alignment/bedgraph/${histControl}_bowtie2.filtered.sorted.bedgraph norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_control_${peakNature}.peaks
# bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_top0.01_${peakNature}.peaks
# done; done
# fi;done