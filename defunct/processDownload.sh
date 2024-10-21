#!/bin/bash
projPath="/mnt/c/users/mjmka/Documents/Bioinformatics" 
seacr="/mnt/c/users/mjmka/Documents/Bioinformatics/tools/SEACR/SEACR_1.3.sh" 
# mkdir $projPath/R_beds/GSE125540
# mkdir $projPath/peakCalling/SEACR/GSE125540
histName="H3K27ac"

arr=($projPath/data/GSE125540/*)
for ((i=0;i<=(${#arr[@]}-1);i++));do
# echo ${arr[i]}
# if [[ ${arr[i]} == GSM4552810_Cech-991-1051_1_S1_R1_001_trimmed_nodups.sorted.bw || "${arr[i]}" == GSM4552810_Cech-991-1051_1_S2_R1_001_trimmed_nodups.sorted.bw ]]
# then histName="INPUT" ;else histName="H3K27me3"; fi
# bigWigToBedGraph ${arr[i]} $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.bedgraph

# echo "filtering bedgraph"
# grep -E "chr1$(printf '\t')|chr2$(printf '\t')|chr3$(printf '\t')|chr4$(printf '\t')|chr5$(printf '\t')|chr6$(printf '\t')\
# |chr7$(printf '\t')|chr8$(printf '\t')|chr9$(printf '\t')|chr10$(printf '\t')|chr11$(printf '\t')|chr12$(printf '\t')|chr13$(printf '\t')\
# |chr14$(printf '\t')|chr15$(printf '\t')|chr16$(printf '\t')|chr17$(printf '\t')|chr18$(printf '\t')|chr19$(printf '\t')|chr20$(printf '\t')\
# |chr21$(printf '\t')|chr22$(printf '\t')|chrX$(printf '\t')|chrY$(printf '\t')|chrM$(printf '\t')" $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.bedgraph >$projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.refiltered.bedgraph


 # echo "sorting bedgraph" ${i}
 # sort -k1,1 -k2,2n $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.refiltered.bedgraph >$projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.sorted.bedgraph
# echo "filtering bedgraph" ${i}
# # Rscript $projPath/scripts/removeBlacklist.R data/ncbi/hg38-blacklist.v2.bed R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.sorted.bedgraph R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.badEOL.bedgraph
# bedtools intersect -v -a $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.sorted.bedgraph -b $projPath/data/ncbi/hg38-blacklist.v2.bed > $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.bedgraph
# # cat -A $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.badEOL.bedgraph| tr [:blank:] \\t | sed "s/.\{0,1\}$//; /^$/d" >$projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.bedgraph
echo "calling peaks" ${i}
# if [[ ${histName} == "H3K27me3" ]];then
#for peakNature in stringent relaxed; do
peakNature="relaxed"
bash $seacr $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/GSE125540/IPSC_${histName}_GSE125540_${i}_seacr_top0.01_${peakNature}.peaks
done
# fi


done

# rm $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.bedgraph
# rm $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.sorted.bedgraph
# rm $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.badEOL.bedgraph
# rm $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.bedgraph
# rm $projPath $projPath/R_beds/GSE125540/IPSC_${histName}_GSE125540_${i}.filtered.bedgraph

# control="IPSC_H3K27me3_GSE125540_4.filtered.bedgraph"

# sample=("IPSC_H3K27me3_GSE125540_0 IPSC_H3K27me3_GSE125540_1 IPSC_H3K27me3_GSE125540_2")
# for peakNature in stringent relaxed; do
# for histName in ${sample[@]};do
# echo "calling peaks" $histName
# bash $seacr $projPath/R_beds/GSE125540/${histName}.filtered.bedgraph 0.01 norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_top0.01_${peakNature}.peaks
# bash $seacr $projPath/R_beds/GSE125540/${histName}.filtered.bedgraph $projPath/R_beds/GSE125540/${control} norm $peakNature $projPath/peakCalling/SEACR/${histName}_seacr_control_${peakNature}.peaks
# done; done