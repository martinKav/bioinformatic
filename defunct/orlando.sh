#!/bin/bash
projPath="/mnt/c/users/mjmka/Documents/Bioinformatics" 
seacr="/mnt/c/users/mjmka/Documents/Bioinformatics/tools/SEACR/SEACR_1.3.sh" 
# mkdir $projPath/R_beds/IPSC_yH2AX_orlando
# mkdir $projPath/peakCalling/SEACR/IPSC_yH2AX_orlando


# bigWigToBedGraph $projPath/data/orlando/GSM5016943_GFP-Input_S1.ucsc.bigWig $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_INPUT_orlando.bedgraph
# seacr="/mnt/c/users/mjmka/Documents/Bioinformatics/tools/SEACR/SEACR_1.3.sh" 

 # echo "filtering bedgraph"
# # grep -E "chr1$(printf '\t')|chr2$(printf '\t')|chr3$(printf '\t')|chr4$(printf '\t')|chr5$(printf '\t')|chr6$(printf '\t')\
# # |chr7$(printf '\t')|chr8$(printf '\t')|chr9$(printf '\t')|chr10$(printf '\t')|chr11$(printf '\t')|chr12$(printf '\t')|chr13$(printf '\t')\
# # |chr14$(printf '\t')|chr15$(printf '\t')|chr16$(printf '\t')|chr17$(printf '\t')|chr18$(printf '\t')|chr19$(printf '\t')|chr20$(printf '\t')\
# # |chr21$(printf '\t')|chr22$(printf '\t')|chrX$(printf '\t')|chrY$(printf '\t')|chrM$(printf '\t')" $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_INPUT_orlando.bedgraph >$projPath/R_beds/IPSC_yH2AX_orlando/IPSC_INPUT_orlando.refiltered.bedgraph


 # # # echo "sorting bedgraph" ${i}
 # # # sort -k1,1 -k2,2n $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_yH2AX_orlando.refiltered.bedgraph >$projPath/R_beds/IPSC_yH2AX_orlando/IPSC_yH2AX_orlando.sorted.bedgraph
# # echo "filtering bedgraph" ${i}
# # # Rscript $projPath/scripts/removeBlacklist.R data/ncbi/hg38-blacklist.v2.bed R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.sorted.bedgraph R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.filtered.badEOL.bedgraph
# # bedtools intersect -v -a $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_INPUT_orlando.refiltered.bedgraph -b $projPath/data/ncbi/hg38-blacklist.v2.bed > $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_INPUT_orlando.filtered.bedgraph
# # cat -A $projPath/R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.filtered.badEOL.bedgraph| tr [:blank:] \\t | sed "s/.\{0,1\}$//; /^$/d" >$projPath/R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.filtered.bedgraph
# echo "calling peaks" ${i}

# for peakNature in stringent relaxed; do

# bash $seacr $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_yH2AX_orlando.filtered.bedgraph $projPath/R_beds/IPSC_yH2AX_orlando/IPSC_INPUT_orlando.filtered.bedgraph norm $peakNature $projPath/peakCalling/SEACR/IPSC_yH2AX_orlando/IPSC_yH2AX_orlando_seacr_control_${peakNature}.peaks
# done
# done

# # rm $projPath/R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.bedgraph
# # rm $projPath/R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.sorted.bedgraph
# # rm $projPath/R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.filtered.badEOL.bedgraph
# # rm $projPath/R_beds/GSE125540/IPSC_H3K27ac_GSE125540_${i}.filtered.bedgraph
 gunzip peakCalling/choi2020/* $projPath/peakCalling/choi2020/GSM4743060_IPSC_H3K27me3_merged_sorted_noDup_peaks.narrowPeak

findPeaks $projPath/R_beds/IPSC_yH2AX_orlando/yH2AX -i $projPath/R_beds/IPSC_yH2AX_orlando/INPUT -style histone -o $projPath/peakCalling/Homer/IPSC_yH2AX_orlando
$projPath/R_beds/IPSC_yH2AX_orlando/yH2AX/IPSC_yH2AX_orlando.refiltered.bed

histName="yH2AX"
sed -e s/\\t/\\t*\\t/3 $projPath/R_beds/IPSC_yH2AX_orlando/${histName}/IPSC_${histName}_orlando.refiltered.bed >$projPath/R_beds/IPSC_yH2AX_orlando/${histName}/IPSC_${histName}_orlando.reformatHomer.bed
makeTagDirectory $projPath/R_beds/IPSC_yH2AX_orlando/${histName} -force5th $projPath/R_beds/IPSC_yH2AX_orlando/${histName}/IPSC_${histName}_orlando.reformatHomer.bed -format bed
