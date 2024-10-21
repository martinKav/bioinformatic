Please Do:
Install the required programs in linux command line and R packages in Rstudio.
Set up your working directory.
Rename the basespace illumina sequencing data folders to your sample names found in histList below. eg this may involve changing the yH2Ax, 2-only, and iPSC prefix. The names of the fastq files themselves don't matter, just the folders.
histList=("IPSC_H3K27ac" "D7_H3K27ac" "IPSC_H3K27me3" "IPSC_yH2AX" "IPSC_2O"  "D7_H3K27me3" "D7_yH2AX" "D7_2O" "Undetermined") 
Move your terminal to your working directory eg using "cd $projPath"
Change your "projPath"

Activate your conda enviroment
Run the CUTandTagPipeline. If you do encounter an error, it will probably be at the stage where it calls the "removeBlackList.R" due to the awkwardness of calling R from the commandLine, so just a heads up.

Working Directory Structure:
I would recommend setting up your folders inside your working directory a manner like:
alignment	this holds the alignments you produce but should be produced automatically by CUTandTagPipeline.sh
data	this holds your reads grabbed straight from illumina but with renamed folders, and also any of the ncbi data like chromosome sizes, index files for bowtie2, and also genehancer or encode data
fastq	this holds your trimmed fastq files
fastqFileQC	this holds the fastqc file reports - it probably should be in the fastq folder but this was just how it was setup originally
MK3	this folder contains the renamed fastq files 
peakCalling	this holds your called peaks
plots	this holds your 
scripts	this is where you keep your script files for organisation, and have a subscript folder for smaller automatic scripts
tools	this is where you install PICARD and SEACR

Script Explanations:
CUTandTagPipeline.sh
This script is our pipeline for our CUT&Tag fastq files. 
Prior to its execution, you need to open it and modify the working directory variable named "projPath".
You also need to double check the variables of "picardCMD", "seacr", "ref" and chromSize to ensure they hold the filepaths for picard/build/libs/picard.jar, SEACR/SEACR_1.3.sh, you bowtie2 index files for hg38 and your hg38 chromosome sizes.
Depending on your file path structure, they may be different so checking is important.
Ensure your blacklist bed file is found at data/ncbi/hg38-blacklist.v2.bed 

processEncode.sh
This script converts the encode bigBeds stored in data/encode/ into beds, renames them to show that MACS2 called their peaks and puts them into the peakCalling folder, preserving their path structure.
It is highly recommended that the encode data is structured to separate the samples by their cell type, histone type and sample name - eg H1/H1_H3K27ac/blablabla

removeBlackList.R
This is a script designed to be called by the command line to remove the blacklisted regions from bed files. Passing filenames to Rscripts via the command line is awkward in terms of relative paths so the "here" library is used.

scaleBedgraphs.R 
This is a script that scales the bedraph files so that their 0.999 quantile score is equal to 1. 

zhengSequencingQC.R
This is a script that generates plots to illustrate the quality of our sample and is adapted from the Ye Zheng C&T Processing Tutorial

PeakPreparation.R
This script calls the filter_peakFileName and annotateAndEnrich functions.

filter_peakFileName.R 
This is a source file for a non-returning function that filters blacklisted regions and alt contigs from peak files, and also converts non-seacr peaks files like those from ENCODE (MACS2) to a seacr-like format.

annotateAndEnrich.R 
This is a source file to annotate a GRanges object to the KnownGene UCSC database, with optional GREAT enrichment (but not recommmended to do GREAT enrichment here as it is better done during the intersection step).

findInDirList.R 
This is a source file for a function that returns a vector of file names in a directory that include a string. It was used to easily grab filenames of peaks files etc

GenehancerPreparation.R 
This script prepares the genehancer dataset by removing blacklisted regions and then annotating our samples to it, producing plots and adding interactions.

EnhancerAtlasPreparation.R 
This script similarly prepares the enhancer atlas dataset by removing blacklisted regions and then annotating our samples to it, producing plots and adding interactions. However, it didn't turn out as well as genehancer so the data was not used.

getOverlapsAndEnrich.R 
This file uses the functions defined in other files to enrich overlaps between samples and also to generate venn diagrams and upset plots (although the upset plots are not good) showing intersections.

buildXwiseListFromVec.R 
This is a source file for a function that returns a list of combinations of components of a string vector. A number of elements can be inputted to choose how many samples to compare at once.

getGRListIntersections.R 
This is a source file to calculate intersections of samples in a GRangesList. 
First it creates a matrix with "in" samples and "out" samples. 
Then it finds what's common in the "in" samples and then keeps what's not found in the "out" samples.

enrichGRange.R 
This is a source file for a function that peforms GREAT, KEGG and enrichGO enrichment on our intersected data and also outputs a condensed form that I shared with Jess.

smallTools.R
This is a collection of one-line functions that I used in several scripts: filterByScoreGR, filterByTotalSignal, mcolsRemove

Other scripts available on request:
I had other scripts to process GEO data like IPSC H3K27ac and H3K27me3 data but didn't end up using it.
I had a script to prepare the Orlando yH2AX data but it's only one sample so I didn't bother including it. It's just bigWigToBedGraph, filtering, sorting, then peak calling with HOMER

Required Command Line programs
cutadapt
fastqc
bowtie2
samtools
picard
bedtools
bedGraphToBigWig
bigWigToBedGraph
bigBedToBed
SEACR
HOMER 
conda is required to install several of these programs so become familiar with it and always install in a virtual environment

Required R packages - hopefully I didn't forget anything
library(stringr)
library(GenomicRanges)
library(here)
library(readr)
library(dplyr)
library(ggplot2)
library(viridis) 
library(ggpubr)
library(plyr)
library(venn)
library(UpSetR)
library(ggplot2)
library(ggpolypath)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(rentrez)
library(rGREAT)
library(R.utils)
library(gtools)

Online file sources:
ENCODE - https://www.encodeproject.org/  go to the stem cell differentiation section and download the replicated peaks files for the H9, H1, Nephron progrenitor and Kidney organoid histones of interest.
Orlando yH2AX data - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5016940 for yH2AX, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5016943 for INPUT, 
hg38 chromosome sizes - https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
genehancer data - https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=2315810080_bdtyDl7auMu8X8mGW2xau99z9V6i&clade=mammal&org=Human&db=hg38&hgta_group=regulation&hgta_track=geneHancer&hgta_table=geneHancerRegElementsDoubleElite&hgta_regionType=range&position=chr1%3A1-248%2C956%2C422&hgta_outputType=primaryTable&hgta_outFileName=
^ you will have to download it chromosome by chromosome and paste them together(notepad++ worked fine but R is probably better)
enhancerAtlas - http://enhanceratlas.org/downloadv2.php
SEdb - https://bio.liclab.net/sedb/  - this site is really laggy by the way just a heads up
bowtie2indexes for hg38 noalt- https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
hg38 blacklist - https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
