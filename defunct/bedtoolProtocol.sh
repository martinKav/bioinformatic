#!/bin/bash
projPath="/mnt/c/users/mjmka/Documents/Bioinformatics" 
chromSize="/mnt/c/users/mjmka/Documents/Bioinformatics/data/ncbi/hg38_chromSizes.fix.txt"
sample=("IPSC_H3K27ac" "IPSC_H3K27me3" "IPSC_yH2AX" "IPSC_2O" "D7_H3K27ac" "D7_H3K27me3" "D7_yH2AX" "D7_2O" "Undetermined")
plusminusList=("10000" "15000")
# "2000" "5000" "500" "0")
typeList=("Enhancer" "PromoterEnhancer")
for type in ${typeList[@]};do
for plusminus in ${plusminusList[@]};do
###### expands at regions plus +- the start site to look at them - we will also expand it in the start coord direction because that's what we observed when graphing
bedtools slop -r 0 -l ${plusminus} -i $projPath/data/Genehancer/tsvs/genehancer_${type}_only.tsv -g $chromSize >$projPath/data/Genehancer/window_beds/genehancer_${type}_plusminus${plusminus}bp_start_extended.bed
###### splits up regions into 5bp windows
bedtools makewindows -b $projPath/data/Genehancer/window_beds/genehancer_${type}_plusminus${plusminus}bp_start_extended.bed -w 5 -i srcwinnum | sort -k1,1 -k2,2n | tr "_" "\t" >$projPath/data/Genehancer/window_beds/genehancer_${type}_plusminus${plusminus}bp_start_extended.5bp.windows.bed
rm $projPath/data/Genehancer/window_beds/genehancer_${type}_plusminus${plusminus}bp_start_extended.bed
for histName in ${sample[@]}; do
bedtools map -a $projPath/data/Genehancer/window_beds/genehancer_${type}_plusminus${plusminus}bp_start_extended.5bp.windows.bed -b $projPath/alignment/bedgraph/${histName}_bowtie2.filtered.sorted.bedgraph -c 4 -o mean -null 0 >$projPath/data/Genehancer/plots/${histName}_genehancer.window.coverage.bedgraph 
sort -t$'\t' -k5,5n $projPath/data/Genehancer/plots/${histName}_genehancer.window.coverage.bedgraph | bedtools groupby -i - -g 5 -c 6 -o sum > $projPath/data/Genehancer/plots/${histName}_genehancer_${type}_plusminus${plusminus}bp_start_extended.window.counts.txt
rm $projPath/data/Genehancer/plots/${histName}_genehancer.window.coverage.bedgraph
done
rm $projPath/data/Genehancer/window_beds/genehancer_${type}_plusminus${plusminus}bp_start_extended.5bp.windows.bed
done
done