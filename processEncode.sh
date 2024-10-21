#!/bin/bash
projPath="/mnt/c/users/mjmka/Documents/Bioinformatics" #change this to your working directory

####### this loop looks for files in a directory called data/encode/ and converts them to bed from bigBed and also renames them to indicate what peak caller was used
####### it preserves the file path structure so having your data structured prior to use is important (eg H1/H1_H3K27ac/ENSMEIFJIEJJI.bigBed)
# get file names
arr=($projPath/data/encode/*)
for ((i=0;i<=${#arr[@]};i++));do
	if [[ -d ${arr[i]} ]]; then
		arr2=(${arr[i]}/*)
		for((j=0;j<=${#arr2[@]};j++));do
			if [[ -d ${arr2[j]} ]]; then
				arr3=(${arr2[j]}/*)
				for((k=0;k<=${#arr3[@]};k++));do
					if [[ -d ${arr3[k]} ]]; then 
						arr4=(${arr3[k]}/*)
						for((l=0;l<=${#arr4[@]};l++));do
							if [[ -d ${arr4[l]} ]];then 
								echo "there's a directory 5 levels deep?" ${arr4[l]}
							else
							
							mkdir -p $projPath/peakCalling/encode/${arr[i]##*/}/${arr2[j]##*/}/${arr3[k]##*/}
							temp=${arr4[l]##*/}
							bigBedToBed ${arr4[l]} $projPath/peakCalling/encode/${arr[i]##*/}/${arr2[j]##*/}/${arr3[k]##*/}/${temp%%.*}_MACS2_peaks.bed
							fi
						done
					else
						if [[ -d $projPath/peakCalling/encode/${arr[i]##*/}/${arr2[j]##*/} ]];then 
							temp=${arr3[k]##*/}
							bigBedToBed ${arr3[k]} $projPath/peakCalling/encode/${arr[i]##*/}/${arr2[j]##*/}/${temp%%.*}_MACS2_peaks.bed
						else 
							mkdir -p $projPath/peakCalling/encode/${arr[i]##*/}/${arr2[j]##*/} 
							temp=${arr3[k]##*/}
							bigBedToBed ${arr3[k]} $projPath/peakCalling/encode/${arr[i]##*/}/${arr2[j]##*/}/${temp%%.*}_MACS2_peaks.bed
						fi 
					fi
				done
			else
				echo ${arr2[j]}
			fi
		done
	else
		echo ${arr[i]}
	fi
done