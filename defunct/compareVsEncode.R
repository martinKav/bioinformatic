# Title: Compare our data to ENCODE
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
# ***************************************************
# I want to find out how many of the peaks are shared between each sample

encode_GRlist <- readRDS(paste0(projPath, "/R_beds/encode/ENCODE_GRanges.GRangesList"))

mxGap=1000
scoreThresholdList= c("top10percent", "top25percent")

for (scoreThreshold in scoreThresholdList)
{
  ourGR_list <- readRDS(file=paste0(projPath, "/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_", scoreThreshold ,"_controlgone_FIXED.peaks.GRangesList"))  
  for (wantedHist in c("H3K27ac", "H3K27me3"))
  {
    # init the matrix
    sampleList <- c(names(ourGR_list[grepl(wantedHist, names(ourGR_list), fixed = TRUE)]), names(encode_GRlist[grepl(wantedHist, names(encode_GRlist), fixed = TRUE)]))
    compDF <- matrix(nrow = length(sampleList), ncol = length(sampleList))
    colnames(compDF) <- sampleList
    rownames(compDF) <- sampleList
    newSampleList <-sampleList
    for (histName in newSampleList)
    {
      if(length(ourGR_list[[histName]]))
      {
        histTest <- ourGR_list[[histName]]
        for (histCompare in newSampleList)
          {
            if(length(ourGR_list[[histCompare]]))
            {
              histControl <- ourGR_list[[histCompare]] 
              compDF[histName, histCompare] <- length(subsetByOverlaps(histTest, histControl, maxgap = mxGap))/length(histTest)
              compDF[histCompare, histName] <- length(subsetByOverlaps(histControl, histTest, maxgap = mxGap))/length(histControl)
            }
            else if(length(encode_GRlist[[histCompare]]))
            {
              histControl <- encode_GRlist[[histCompare]]
              testInCompare=0
              compareInTest=0
              for (subHistCompare in names(histControl))
              {
                # length of the subset/length of the test 
                histSubControl <- histControl[[subHistCompare]]
                testInCompare = testInCompare + (length(subsetByOverlaps(histTest, histSubControl, maxgap = mxGap))/length(histTest))
                compareInTest = compareInTest + (length(subsetByOverlaps(histSubControl, histTest, maxgap = mxGap))/length(histSubControl))
              }
              compDF[histName, histCompare] <- testInCompare/length(histControl)
              compDF[histCompare, histName] <- compareInTest/length(histControl)
            }
            else {print("errror finding histControl")}
          }
      }
      else if(length(encode_GRlist[[histName]]))
      {#if its not ours, get it from encode
        histTest <- encode_GRlist[[histName]]
        
          for (histCompare in newSampleList)
          {
            if(length(ourGR_list[[histCompare]]))
            {
              histControl <- ourGR_list[[histCompare]] 
              testInCompare=0
              compareInTest=0
              for (subHistName in names(histTest))
              {
                histSubTest <- histTest[[subHistName]]
                testInCompare = testInCompare + (length(subsetByOverlaps(histSubTest, histControl, maxgap = mxGap))/length(histSubTest))
                compareInTest = compareInTest + (length(subsetByOverlaps(histControl, histSubTest, maxgap = mxGap))/length(histControl))
                
              }
              compDF[histName, histCompare] <- testInCompare/length(histTest)
              compDF[histCompare, histName] <- compareInTest/length(histTest)
            }
            else if(length(encode_GRlist[[histCompare]]))
            {
              histControl <- encode_GRlist[[histCompare]]
              testInCompare=0
              compareInTest=0
              for (subHistName in names(histTest))
              {
                histSubTest <- histTest[[subHistName]]
                for (subHistCompare in names(histControl))
                {
                  histSubControl <- histControl[[subHistCompare]]
                  testInCompare = testInCompare + (length(subsetByOverlaps(histSubTest, histSubControl, maxgap = mxGap))/length(histSubTest))
                  compareInTest = compareInTest + (length(subsetByOverlaps(histSubControl, histSubTest, maxgap = mxGap))/length(histSubControl))
                }
                
              }
              compDF[histName, histCompare] <- testInCompare/(length(histTest)*length(histControl))
              compDF[histCompare, histName] <- compareInTest/(length(histTest)*length(histControl))
              
            }
            else {print("errror finding histControl")}
          }
      }
      else {print("errror finding histTest")}
      
      #remove histName from sampleList
      
      newSampleList <- sampleList[sampleList != histName]
    }
    write.table(compDF, file=paste0(projPath, "/R_beds/bulkData/", wantedHist, "_", scoreThreshold, "_crossComparison_", mxGap,"BPgap.txt"))
  }
}
wantedHist="H3K27ac"
scoreThreshold="top25percent"
compDF <- read.table(file=paste0(projPath, "/R_beds/bulkData/", wantedHist, "_", scoreThreshold, "_crossComparison_", mxGap,"BPgap.txt"))
comp_pca <- prcomp(compDF, center=TRUE, scale=TRUE)
summary(comp_pca) 
comp_pcaDF <- data.frame(comp_pca$rotation)
compDF$names <- row.names(compDF)
compDF[grepl("H1", names(compDF), fixed = TRUE),"colors"]  <- "red"   
compDF[grepl("H9", names(compDF), fixed = TRUE),"colors"]  <- "green"   
compDF[grepl("Nephron", names(compDF), fixed = TRUE),"colors"]  <- "pink"   
compDF[grepl("NPC", names(compDF), fixed = TRUE),"colors"]  <- "pink"
compDF[grepl("IPSC", names(compDF), fixed = TRUE),"colors"]  <- "black"
compDF[grepl("D7", names(compDF), fixed = TRUE),"colors"]  <- "blue"   
comp_pcaDF

library(colorBlindness)
library(viridis)
p <-prcomp(overlapMatrix, center=TRUE, scale=TRUE)
p
summary(p)
pDF$names <- row.names(p)
length(pDF)
pDF$colors <- rainbow(length(pDF))
pDF[,1:54]
colnames(pDF)
library(cluster)
library(factoextra)
fviz_nbclust(p$rotation, kmeans, method = "silhouette")
pDF$cluster <- kmeans(p$rotation, centers= 3, iter.max=100,nstart=10, algorithm = "MacQueen")$cluster
pDF$cluster


library(ggfortify)
ggplot(data=p, mapping = aes(x=PC1, y= PC2))+
  geom_point(colour=pDF$colors, size = 2)
biplot_p <- biplot(p)

library(threejs)
z <- p$rotation[,"PC3"]
x <- p$rotation[,"PC1"]
y <- p$rotation[,"PC2"]
scatterplot3js(x,y,z, color=pDF$colors)

library(plotly)
plot_ly(pDF, x = pDF[,"PC1"], y = pDF[,"PC2"], 
        text = paste("Clarity: ", row.names(pDF)),
        mode = "markers", color = pDF$colors, size = pDF[,"PC3"])

qDF <- data.frame(q$rotation)
qDF$names <- row.names(qDF)
qDF$colors <- c("red", "orange", "blue", "green", "pink", "purple", "cyan", "magenta", "brown")


png(file=paste0(projPath, "/plots/PCA1.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
library(ggfortify) 
ggplot(data=q,mapping=aes(x=PC1, y=PC2))+
  geom_point(colour=qDF$colors, size=6)
dev.off()

png(file=paste0(projPath, "/plots/PCA2.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
biplot_comp <- biplot(q)  
dev.off()

png(file=paste0(projPath, "/plots/PCA3.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
library(threejs)
z <- q$rotation[,"PC3"]
x <- q$rotation[,"PC1"]
y <- q$rotation[,"PC2"]
scatterplot3js(x,y,z, color=qDF$colors)
dev.off()

png(file=paste0(projPath, "/plots/PCA4.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
library(plotly)
plot_ly(qDF, x = qDF[,"PC1"], y = qDF[,"PC2"], 
        text = paste("Clarity: ", row.names(qDF)),
        mode = "markers", color = qDF$colors, size = qDF[,"PC3"])
dev.off()


library(ggfortify) 
ggplot(data=q,mapping=aes(x=PC1, y=PC2))+
  geom_point(colour=qDF$colors, size=6)
  
biplot_comp <- biplot(comp_pca)    

library(threejs)
z <- comp_pca$rotation[,"PC3"]
x <- comp_pca$rotation[,"PC1"]
y <- comp_pca$rotation[,"PC2"]
scatterplot3js(x,y,z, color=compDF$colors)

library(plotly)
plot_ly(comp_pcaDF, x = comp_pcaDF[,"PC1"], y = comp_pcaDF[,"PC2"], 
        text = paste("Clarity: ", row.names(comp_pcaDF)),
        mode = "markers", color = compDF$colors, size = comp_pcaDF[,"PC3"])

 # ok so this is clearly a sign that iPSCs are not good equivalents to ESCs
# lets try find some data online

"~/Bioinformatics/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_top10percent_controlgone_FIXED.peaks.GRangesList"
"~/Bioinformatics/R_beds/CUTandTAG/filtered/sample_seacr_relaxed_top25percent_controlgone_FIXED.peaks.GRangesList"
"~/Bioinformatics/R_beds/encode/ENCODE_GRanges.GRangesList"

blackBed <- read.table(file = paste0(projPath, "/R_beds/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
colnames(blackBed) <- c("chrom", "start", "end", "reason")
blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
                          seqnames=blackBed$chrom,
                          reason=blackBed$reason)



# is A in B the same as B in A? no!

# test granges just to confirm overlap behaviour
a <- data.frame(chr="chr1",
                start=c(1,10000,20000),
                end=c(50, 10020, 25000))
a
b <- data.frame(chr = "chr1",
                start=c(10,10010,20100),
                end=c(40,10200,21000))
grA <- GRanges( seqnames=a$chr,
                ranges=IRanges(start=a$start,
                               end=a$end))
grB <- GRanges( seqnames=b$chr,
                  ranges=IRanges(start=b$start,
                                 end=b$end))
subsetByOverlaps(grA, grB)
subsetByOverlaps(grA, grA, maxgap = 100000)

findOverlaps(grA, grB)











c(
  "~/Bioinformatics/R_beds/genehancer/GRangesLists/genehancer_filtered.GRangesList",
  "~/Bioinformatics/R_beds/genehancer/GRangesLists/genehancer_filtered_500plus.GRangesList"
)
c(
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_hg38_ele_filtered_500plus.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_hg38_filtered_500plus.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_hg38_ele_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_hg38_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Renal_hg38_ele_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Renal_hg38_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Kidney_hg38_ele_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Kidney_hg38_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_iPSCs_hg38_ele_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_iPSCs_hg38_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Embryonic_kidney_hg38_ele_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Embryonic_kidney_hg38_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Embryo_hg38_ele_filtered.GRangesList",
  "~/Bioinformatics/R_beds/SEdb/GRangesLists/Human_Embryo_hg38_filtered.GRangesList"
)