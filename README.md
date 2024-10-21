I wrote these scripts to analyse CUT&Tag sequencing data for 3 histone modifications: H2AX S139P, H3K27me3, and H3K27ac. We wanted to compare the localisation of H2AX S139P to enhancers involved in cell fate across different cell types. Enhancer databases were downloaded and compared in Rstudio. 

We were unable to get any intelligible results due to the large amount of off-target transposase activity in our sample. Additionally, our read counts were very low for some of our samples, and the majority of our reads were unable to be demultiplexed.

These scripts are not very performant. I was under a lot of time pressure when doing the bioinformatic analysis because of how late I got my results. If I had to do it again, I would make much heavier use of data.table or write some of the intensive scripts in rcpp. 
