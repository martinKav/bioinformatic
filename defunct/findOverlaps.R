# TItle: Find Overlaps in Cut and Tag peaks
# Author: Martin Kavanagh
# Date: 01/08/24
# ***************************************************
rm(list=ls())
library(GenomicRanges)
library(stringr)
library(readr)
library(ggplot2)
setwd("C:/Users/mjmka/Documents/Bioinformatics")
projPath="C:/Users/mjmka/Documents/Bioinformatics"
# **************************************************
# we'll just start with our data 
# read GRangesLists
# calc overlap
peakNature="relaxed"
statsTableList <- list()
scoreThresholdList= c("top10percent", "top25percent")
for (scoreThreshold in scoreThresholdList)
{
  ourGRlist <- readRDS(file=paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_",scoreThreshold,"_controlgone_Filtered.peaks.GRangesList"))
  statsTable <- data.frame(row.names = names(ourGRlist), col.names = names(ourGRlist))
  for (histName in names(ourGRlist))
  {
    for (histCompare in names(ourGRlist))
    {
      histTest <- unlist(ourGRlist[histName])
      histControl <- unlist(ourGRlist[histCompare])
      overlap <- findOverlaps(histTest, histControl, type="any", select="all", ignore.strand = T)
      statsTable[histName, histCompare] <- length(overlap)/length(histTest)
    }
  }
  statsTableList[[scoreThreshold]] <- statsTable
}
statsTableList

sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX", 
               #"IPSC_2O", "D7_2O", 
               "Undetermined")
      
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX", "Undetermined")
                              
# creating the top 0.25 and top 0.1
peakNatureList=c("relaxed","stringent")
for (peakNature in peakNatureList)
{
  #removing IPSC_2O and D7_2O
  ourGRlist <- readRDS(file=paste0(projPath, "/R_beds/CUTandTAG/sample_top0.01.peaks.",peakNature,".GRangesList"))
  gr_vector <- list()
  for (histName in sampleList)
  {
    histTest <- ourGRlist[[histName]]
    for (histCompare in c("IPSC_2O", "D7_2O"))
    {
      histControl <- ourGRlist[[histCompare]]
      histTest <- subsetByOverlaps(histTest, histControl, type="any", invert=TRUE)
    }
    gr_vector[[histName]] <- histTest
  }
  ourGR_list <- GRangesList(gr_vector)
  # filtering score
  scoreThresholdList= c("top10percent", "top25percent", "top50percent")
  for (scoreThreshold in scoreThresholdList)
  {
    gr_vector <- list()
    for (histName in names(ourGR_list))
    {
      gr <- ourGR_list[[histName]]
      gr
      quantile(score(gr), probs=c(0.9,0.1))[1]
      PROBS=c()
      if(scoreThreshold == "top10percent"){PROBS=c(0.9,0.1)} 
      else if(scoreThreshold=="top25percent"){PROBS=c(0.75,0.25)}
      else{PROBS=c(0.5,0.5)}
      grG <- gr[score(gr) >quantile(score(gr), probs=PROBS[1])]
      gr_vector[[histName]] <- grG 
      
    }
    GRlist_new <- GRangesList(gr_vector)
    saveRDS(GRlist_new, file=paste0(projPath, "/R_beds/CUTandTag/filtered/sample_seacr_", peakNature, "_", scoreThreshold, "_Filtered.peaks.GRangesList"))
  }
}

# find overlaps between them
# get stats about overlap 
# amount of overlaps, percentage overlap, how much overlap on average, fully overlap, how much overlap if 1000bp dist given
#get overlap stats
row_list=(c(rep("IPSC_H3K27me3", 6), rep("IPSC_H3K27ac",6), rep("IPSC_yH2AX", 6), rep("D7_H3K27me3",6), rep("D7_H3K27ac",6), rep("D7_yH2AX",6),rep("Undetermined",6)))
row_list=1:(length(names(ourGR_list))*(length(names(ourGR_list))-1))
col_list = c("NumOverlap", "PctOverlap", "MeanOverlap","NumHighOverlap","PctHighOverlap", "NumPeaksDistant", "PctPeaksDistant")
peakNature="relaxed"
#for (peakNature in peakNatureList)
{
  scoreThresholdList= c("top10percent", "top25percent", "top50percent")
  for (scoreThreshold in scoreThresholdList)
  {
    
    ourGR_list <- readRDS(file=paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_",peakNature,"_",scoreThreshold,"_Filtered.peaks.GRangesList"))
    row_list=1:(length(names(ourGR_list))*(length(names(ourGR_list))-1))
    overlapData <- matrix(nrow = length(row_list), ncol = length(col_list))
    colnames(overlapData) <- col_list 
    rownames(overlapData) <- c(1:nrow(overlapData))
    i=1
      for (histName in names(ourGR_list))
      {
        histTest <- ourGR_list[[histName]]
        for (histControl in names(ourGR_list))
        {
          histCompare <- ourGR_list[[histControl]]
          
            if(histName != histControl)
            {
              overlap <- findOverlaps(histTest, histCompare,type="any", select="all", ignore.strand = T)
              overlapData[i,"NumOverlap"] <- length(overlap)
              overlapData[i,"PctOverlap"] <- length(overlap)/length(histTest)
              avg_list <- vector(mode="numeric")
              for(j in 1:length(overlap))
              {
                avg_list[j] <- max(c(end(histTest[queryHits(overlap[j])]) -start(histCompare[subjectHits(overlap[j])]),
                                     end(histCompare[subjectHits(overlap[j])]) - start(histTest[queryHits(overlap[j])])))
              }
                print(paste0("finished getting average overlap of", histName, " vs ", histControl ))
                overlapData[i,"MeanOverlap"] <- mean(avg_list, na.rm=T)
                
                # high though mmmmmmm, will need to calc width of peak
                # mmmmmmm can get width of query, and can do subset to get 
                # build Grange one by one hahahahaha
                # so start, find min width of histname, then subset with that as minoverlap/50
                # then go through the subset 1 by 1 and 
                # get hits object - compare start vs end of other -> OR startment (if end2 - start1 or end1 - start2)
                highOverLapIndex <- findOverlaps(histTest, histCompare, minoverlap=min(width(histTest), na.rm=T)/2 ,type="any", select="all", ignore.strand = T)
                reduced_list <- vector(mode="numeric")
                k=1
                for(j in 1:length(highOverLapIndex))
                {
                  cond <-(end(histCompare[subjectHits(highOverLapIndex[j])]) -start(histTest[queryHits(highOverLapIndex[j])])) > width(histTest)[queryHits(highOverLapIndex[j])]/2 ||
                    end(histTest[queryHits(highOverLapIndex[j])]) - start(histCompare[subjectHits(highOverLapIndex[j])]) > width(histTest)[queryHits(highOverLapIndex[j])]/2
                  if(!is.na(cond) && cond )
                  {
                    reduced_list[k] <- queryHits(highOverLapIndex[j])
                    k = k+1
                  }
                }
                print(paste0("finished big overlaps of", histName, " vs ", histControl ))
                overlapData[i,"NumHighOverlap"] <- length(reduced_list)
                overlapData[i, "PctHighOverlap"] <- length(reduced_list)/length(histTest)
                # cant do highoverlap simply or low overlap simply
                # low overlap actually fine, just invert subset with a maxgap of 2000
                distNoOverlap <- subsetByOverlaps(histTest, histCompare, maxgap =2000, type="any", ignore.strand = T, invert=TRUE)
                overlapData[i,"NumPeaksDistant"] <- length(distNoOverlap)
                overlapData[i,"PctPeaksDistant"] <- length(distNoOverlap)/length(histTest)
                rownames(overlapData)[i] <- paste0(histName,"_",histControl)
                print(paste0("Finished ", histName, " vs ", histControl))
                i =i+1
              }
          }
      
      }
      
    write.table(overlapData, file=paste0(projPath,"/R_beds/CUTandTAG/filtered/sample_seacr_",peakNature,"_", scoreThreshold,", overlap_stats_Filtered.peaks.txt"))
  }
}
names(ourGRlist)
ourGRlist$IPSC_2O
# subset by overlaps
# gets what regions are common in both samples, and what is different
# get granges object, and apply genome ontology enrichment to it


# genome ontology should be easy thereon