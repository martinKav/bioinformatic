#!/bin/bash
projPath="/mnt/c/users/mjmka/Documents/Bioinformatics"  # change this to the path of your working directory
ref="${projPath}/data/bowtie2Index/hg38/GRCh38_noalt_as" # change this to your bowtie2 index files for hg38
picardCMD="java -jar $projPath/tools/picard/build/libs/picard.jar" # this should point to the picard.jar file
cores=8 # how many cores you want bowtie2 to use
chromSize="${projPath}/data/ncbi/hg38_chromSizes.fix.txt" # this should point to the file with your chromosome sizes for hg38
seacr="${projPath}/tools/SEACR/SEACR_1.3.sh"  # this should point to your seacr script file 
histList=("IPSC_H3K27ac" "D7_H3K27ac" "IPSC_H3K27me3" "IPSC_yH2AX" "IPSC_2O"  "D7_H3K27me3" "D7_yH2AX" "D7_2O" "Undetermined") # this is a list of the names of the histone files you want to process - I renamed them to these to make them easier to work with
 for histName in ${histList[@]} ; do # for loop through each sample name in histList
 #####delete previous files from unsucessful runs (you might not want to do this so it is commented out)
# # # rm -f  $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
# # # rm -f  $projPath/alignment/rmDuplicate/picard_summary/${histName}_picard.rmDup.txt
# # # rm -f  $projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.rmDup.txt
# # # rm -f  $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragmentsCount.bin$binLen.bed
# # # rm -f  $projPath/alignment/sam/${histName}_bowtie2.sam
# # # rm -f  $projPath/alignment/sam/${histName}_bowtie2.sorted.sam
# # # rm -f  $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam
# # # rm -f  $projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam
# # # rm -f  $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam
# # # rm -f  $projPath/alignment/bed/${histName}_bowtie2.rmDup.bed
# # # rm -f  $projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed
# # # rm -f  $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed
# # # rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph
# # # rm -f  $projPath/fastq/${histName}_R1.QC.fastq.gz
# # # rm -f  $projPath/fastq/${histName}_R2.QC.fastq.gz

#### make folders - hope I have them all
mkdir $projPath/alignment
mkdir $projPath/alignment/bam
mkdir $projPath/alignment/sam
mkdir $projPath/alignment/bed
mkdir $projPath/alignment/bedgraph
mkdir $projPath/alignment/rmDuplicate
mkdir $projPath/fastq
mkdir $projPath/fastqFileQC
mkdir ${projPath}/alignment/sam/bowtie2_summary
mkdir $projPath/alignment/sam/fragmentLen
mkdir $projPath/alignment/rmDuplicate/picard_summary/
mkdir $projPath/MK3/
mkdir $projPath/MK3/${histName}/
mkdir $projPath/fastqFileQC/${histName}
mkdir $projPath/fastqFileQC/${histName}_postcut
mkdir $projPath/fastq/cut_results/
mkdir $projPath/alignment/bigWig

##### simplify file names for the fasta files
cp $projPath/data/${histName}/*_R1_001.fastq.gz $projPath/MK3/${histName}/${histName}_R1.fastq.gz
cp $projPath/data/${histName}/*_R2_001.fastq.gz $projPath/MK3/${histName}/${histName}_R2.fastq.gz

echo "cutadapting" ${histName} # the cutadapt section of the analysis honestly doesnt change that much but there are a few 100 or so bad reads
cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -j 8 -m 12 -q 25,25 --nextseq-trim 20 -o $projPath/fastq/${histName}_R1.QC.fastq.gz -p $projPath/fastq/${histName}_R2.QC.fastq.gz $projPath/MK3/${histName}/${histName}_R1.fastq.gz $projPath/MK3/${histName}/${histName}_R2.fastq.gz >${projPath}/fastq/cut_results/${histName}_cutadapt_summary.txt

echo "fastq examining cut" ${histName} # simple fastqc after cutadapt
fastqc -o ${projPath}/fastqFileQC/${histName}_postcut -f fastq ${projPath}/fastq/${histName}_R1.QC.fastq.gz
fastqc -o ${projPath}/fastqFileQC/${histName}_postcut -f fastq ${projPath}/fastq/${histName}_R2.QC.fastq.gz

# # # # # # # ######not using bwa but I had tried at one strage
# # # # # # # $projPath/tools/bwa/bwa mem -t 8 -M \
# # # # # # # -R "@RG\tID:${histName}\tPL:ILLUMINA\tSM:${histName}"\
# # # # # # # -P $projPath/bwaIndex/GCF_000001405.40_GRCh38.p14_genomic.fna  \
# # # # # # # $projPath/fastq/${histName}_R1.QC.fastq.gz $projPath/fastq/${histName}_R2.QC.fastq.gz \
# # # # # # # 2> $projPath/alignment/sam/bwa_log/${histName}_bwa.err \
# # # # # # # > $projPath/alignment/sam/${histName}_bwa.sam

echo "aligning with bowtie2" $histName # bowtie2 is our aligner of choice 
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 \
-p ${cores} --rg-id ${histName} --rg SM:${histName} --rg PL:ILLUMINA --rg LB:MK3 --rg PU:MK3_${histName} -x ${ref} -1 \
${projPath}/fastq/${histName}_R1.QC.fastq.gz -2 ${projPath}/fastq/${histName}_R2.QC.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt

echo "sorting sam" $histName ## picard is needed to remove optical duplicates - real duplicates can be removed if needed but it is not important in C&T because duplicate generation is intentional
$picardCMD SortSam -I $projPath/alignment/sam/${histName}_bowtie2.sam -O $projPath/alignment/sam/${histName}_bowtie2.sorted.sam --SORT_ORDER coordinate
echo "removing dupes" ${histName}
$picardCMD MarkDuplicates -I $projPath/alignment/sam/${histName}_bowtie2.sorted.sam -O $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam \
	--REMOVE_SEQUENCING_DUPLICATES true -M $projPath/alignment/rmDuplicate/picard_summary/${histName}_picard.rmDup.txt

echo "grabbing mapped size of dup-removed" $histName # this is for QC analysis to see how many non duplicated fragments were mapped successfully
samtools view -F 0x04 $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort |uniq -c |awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.rmDup.txt
echo "converting to bam" $histName
samtools view -h -bS -F 0x04 $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam >$projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam

##### because the files are so big the intermediate files need to be removed
 echo "removing intermiediate sams" ${histName} 
 rm -f  $projPath/alignment/sam/${histName}_bowtie2.sam
 rm -f  $projPath/alignment/sam/${histName}_bowtie2.sorted.sam
 rm -f  $projPath/alignment/rmDuplicate/${histName}_bowtie2.sorted.rmDup.sam

echo "sorting bams" ${histName} #for some weird reason the bams need to be re-sorted 
samtools sort -n $projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam -o $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam
echo "bam to bed" ${histName}
bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam -bedpe >$projPath/alignment/bed/${histName}_bowtie2.rmDup.bed
echo "cleaning bed" ${histName} # removes fragments with ends on different chromosomes
awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.rmDup.bed >$projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed
echo "extracting frag columns" ${histName} # grabs wanted columns and sorts them
cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed
echo "removing any leftover alt contigs" ${histName} # even though our bowtie2 index should have no alt contigs for hg38, I found that sometimes they slip through for some reason so I just refiltered them 
grep -E "chr1$(printf '\t')|chr2$(printf '\t')|chr3$(printf '\t')|chr4$(printf '\t')|chr5$(printf '\t')|chr6$(printf '\t')\
|chr7$(printf '\t')|chr8$(printf '\t')|chr9$(printf '\t')|chr10$(printf '\t')|chr11$(printf '\t')|chr12$(printf '\t')|chr13$(printf '\t')\
|chr14$(printf '\t')|chr15$(printf '\t')|chr16$(printf '\t')|chr17$(printf '\t')|chr18$(printf '\t')|chr19$(printf '\t')|chr20$(printf '\t')\
|chr21$(printf '\t')|chr22$(printf '\t')|chrX$(printf '\t')|chrY$(printf '\t')|chrM$(printf '\t')" $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed >$projPath/alignment/bed/${histName}_bowtie2.refiltered.bed
# # # # # # # echo "assessing reproductibility" # we would do this if we had replicates but we dont
# # # # # # # binLen=500
# # # # # # # awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.refiltered.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.rmDup.fragmentsCount.bin$binLen.bed

##### because the files are so big the intermediate files need to be removed, but the refiltered bed files and the rmDup mapped bam are kept because they are used in the next step and also in QC
echo "removing intermiediate bams and beds" ${histName}
rm -f $projPath/alignment/bam/${histName}_bowtie2.rmDup.sorted.mapped.bam
rm -f $projPath/alignment/bed/${histName}_bowtie2.rmDup.bed
rm -f  $projPath/alignment/bed/${histName}_bowtie2.rmDup.clean.bed
rm -f  $projPath/alignment/bed/${histName}_bowtie2.rmDup.fragments.bed

echo "getting scale factor" ${histName} # this is a scale factor for scaling our bedraph files and it depends on the library size calculated via the bam file size
TmpScale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c $projPath/alignment/bam/${histName}_bowtie2.rmDup.mapped.bam)")
echo "bed to bedgraph" ${histName}  # produces a bedgraph for visualisation
bedtools genomecov -scale $TmpScale -bg -scale $TmpScale -i $projPath/alignment/bed/${histName}_bowtie2.refiltered.bed -g $chromSize >$projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph
echo "sorting bedgraph" ${histName} # sorts the bedgraph so it can be filtered
sort -k1,1 -k2,2n $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph >$projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.sorted.bedgraph
echo "filtering bedgraph" ${histName} # this an Rscript to filter the blacklisted regions from the bedraph for visualisation - the Rscript is really awkward so it needs to have relative filenames for arguments. I had it setup at one point to accept absolute filenames but then I wiped my data using rm -rf by mistake.
Rscript $projPath/scripts/removeBlacklist.R data/ncbi/hg38-blacklist.v2.bed alignment/bedgraph/${histName}_bowtie2.rmDup.sorted.bedgraph alignment/bedgraph/${histName}_bowtie2.filtered.badEOL.bedgraph
cat -A $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.badEOL.bedgraph| tr [:blank:] \\t | sed "s/.\{0,1\}$//; /^$/d" >$projPath/alignment/bedgraph/${histName}_bowtie2.filtered.bedgraph
echo "sorting and merging bedgraph" ${histName}
sort -k1,1 -k2,2n $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.bedgraph >$projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph

####### after this is just for bigwig generation, honestly its not worth it unless you are strapped for size/want to view using UCSC genome browser so I commented it out - for the visualisation however I instead scaled bedgraph files using scaleBedgraphs.R
# bedtools merge -i $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph -c 4 -d 0 -o max > $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.bedgraph
# LC_COLLATE=C sort -k1,1 -k2,2n $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.bedgraph >$projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.sorted.bedgraph
# echo "converting bedgraph to bigwig" $histName
# bedGraphToBigWig $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.sorted.bedgraph $chromSize $projPath/alignment/bigWig/${histName}.rmDup.bigWig

###### this calls peaks in seacr using the top 0.01 peaks, but with the stringent and relaxed thresholds
for peakNature in stringent relaxed; do
bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_top0.01_${peakNature}.peaks
done

##### removes the intermediate bedgraph files 
echo "removing intermiediate beds etc" ${histName}
rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.fragments.normalized.bedgraph
rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.rmDup.sorted.bedgraph
rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.badEOL.bedgraph 
rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.bedgraph
rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.bedgraph
rm -f  $projPath/alignment/bedgraph/${histName}_bowtie2.merged.max.sorted.bedgraph
 done 

#### this is legacy code from when I was calling peaks using controls, but it is not good because our secondary only controls were so bad, so I switched entirely to calling the top 0.01 peaks
# # # for histControl in IPSC_2O D7_2O; do
# # # if [[ histControl == IPSC_2O ]]; then
# # # for histName in IPSC_H3K27ac IPSC_H3K27me3 IPSC_yH2AX D7_2O Undetermined; do
# # # echo "peakCalling for" histName "with" $histControl
# # # for peakNature in stringent relaxed; do
# # # bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph $projPath/alignment/bedgraph/${histControl}_bowtie2.filtered.sorted.bedgraph norm $peakNature $projPath/peakCalling/SEACR/CUTandTag/${histName}_seacr_control_${peakNature}.peaks
# # # # bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_top0.01_${peakNature}.peaks
# # # done; done
# # # else
# # # for histName in IPSC_2O D7_H3K27ac D7_H3K27me3 D7_yH2AX;do
# # # echo "peakCalling for" histName "with" $histControl
# # # for peakNature in stringent relaxed; do
# # # bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph $projPath/alignment/bedgraph/${histControl}_bowtie2.filtered.sorted.bedgraph norm $peakNature $projPath/peakCalling/SEACR/CUTandTag/${histName}_seacr_control_${peakNature}.peaks
# # # # bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_top0.01_${peakNature}.peaks
# # # done; done
# # # fi;done

####### legacy, used to use a normalisation method but it wouldnt work properly for the underloaded IPSC_H3K27ac sample
# # # echo "executing ChIPseqSpikeInFree script"
# # # rm -f  $projPath/alignment/ChIPseqSpikeInFree 
# # # mkdir $projPath/alignment/ChIPseqSpikeInFree
# # # Rscript $projPath/scripts/ChIPSpikeInFree.R $projPath
# # # rm -f  $projPath/alignment/ChIPseqSpikeInFree/SpikeInFree_results_rawCounts.txt