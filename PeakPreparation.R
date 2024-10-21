# Purpose: to filter our bed files and annotate them.
# Author: Martin Kavanagh
# Date: 01/08/24
# ***************************************************
rm(list=ls())
library(GenomicRanges)
library(stringr)
library(readr)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
setwd(projPath)
source(paste0(projPath, "/scripts/subscripts/filter_peakFileName.R")) # testing the filter but not needed anymore
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
source(paste0(projPath, "/scripts/subscripts/annotateAndEnrich.R"))
# ***************************************************

wantedHists=c("Homer", "encode", "CUTandTag")

fileNames <- vector(mode="character")
fileNames <- findInDirList(path = paste0(projPath, "/peakCalling"), wantedstrings = wantedHists)
fileNames <- fileNames[(!(grepl("stringent",fileNames[file.info(fileNames)$size != 0])))]
fileNames <- fileNames[(!(grepl("2O|Undetermined",fileNames)))]
fileNames
for(fileName in fileNames){
  filter_peakFileName(fileName, remove_non_seacr = TRUE)
}

fileNames <- findInDirList(path = paste0(projPath, "/peakCalling"), wantedstrings = wantedHists)
fileNames <- fileNames[(!(grepl("stringent",fileNames[file.info(fileNames)$size != 0])))]
fileNames <- fileNames[(!(grepl("2O|Undetermined",fileNames)))]
fileNames

for(fileName in fileNames)
{
  annotateAndEnrich(fileName)

}

