# enhancerAtlas Preparation and Investigation
# Author: Martin Kavanagh
# Date: 15/08/2024
# ******************************************
rm(list=ls())
shell("cls")
library(plyr)
library(GenomicRanges)
library(stringr)
library(readr)
library(wrapr)
library(ggplot2)
library(ggsankey)
library(R.utils)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
enhancerSampleList=c("H1", "H9", "IPSC")
colourVec=c("dodgerblue2", "#E31A1C", # red
            "green4",
            "#6A3D9A", # purple
            "#FF7F00", # orange
            "black", "gold1",
            "skyblue2", "#FB9A99", # lt pink
            "palegreen2",
            "cyan", # lt purple
            "#FDBF6F", # lt orange
            "peachpuff3", "khaki2",
            "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
            "darkturquoise", "green1", "yellow4", "yellow3",
            "darkorange4", "brown")
quantList=c(0, 0.5, 0.75, 0.9, 0.99)
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX",
               "H1_H3K27ac","H1_H3K27me3", "H1_H3K4me1", "H9_H3K27ac","H9_H3K27me3", "H9_H3K4me1", "IPSC_ORLANDO_yH2AX")
level_order <- c('Promoter', 'Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', "Promoter (3-4kb)", "Promoter (4-5kb)",
                 "Promoter (5-6kb)", "Promoter (6-7kb)", "3' UTR", "Exonic", "Intronic", "5' UTR", "Downstream (<=300bp)",
                 "Distal Intergenic") 
EH_level_order <- c("At Enhancer", ">50 bp from Enhancer","50-100 bp from Enhancer","200-500 bp from Enhancer","500-1k bp from Enhancer",
                    "1k-2.5k bp from Enhancer","2.5k-5k bp from Enhancer", "5k-10k bp from Enhancer", "<10k bp from Enhancer" )

# data annotated to TSSs with library(TxDb.Hsapiens.UCSC.hg38.knownGene)
blackBed <- read.table(file = paste0(projPath, "/data/ncbi/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
colnames(blackBed) <- c("chrom", "start", "end", "reason")
blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
                          seqnames=blackBed$chrom,
                          reason=blackBed$reason)
enhancerFileNames <- findInDirList(paste0(projPath, "/data/"), "enhancerFiles")
enhancerInteractionFileNames <- findInDirList(paste0(projPath, "/data/"), "InteractionFiles")

fileNames <- findInDirList(paste0(projPath, "/annotationsAndEnrichment/"), c(
  #"encode", "orlando", 
  "CUTandTag"))
fileNames <- fileNames[!(grepl("_2O|Undetermined", fileNames))]
fileNames <- fileNames[!(grepl("enrichment", fileNames))]
fileNames <- fileNames[!(grepl("stringent", fileNames))]
fileNames <- fileNames[!(grepl("enhancerAnnotated", fileNames))]
fileNames

for(enhancerFile in enhancerFileNames)
{
  #trim sides
  enhancerName <- enhancerSampleList[str_detect(enhancerFile, enhancerSampleList)]
  EHbed <- read.table(file=enhancerFile, sep="\t", quote="", header=F)
  EHGR <- GRanges(seqnames=EHbed$V1 ,ranges = IRanges(start=EHbed$V2, end=EHbed$V3),score=EHbed$V4, name=EHbed$V5)
  EHGR <- subsetByOverlaps(EHGR, blacklist.hg38, type = "any", invert = TRUE)
  EHGR
  for(fileName in fileNames)
  {
    #annotate
    for (quant in quantList)
    {
      histName <- sampleList[str_detect(fileName, sampleList)]
      if(!dir.exists(paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName)))
      {
        mkdirs(paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName))
      } 
      annoBed <- read.table(file=fileName, sep="\t", quote="")
      annoBed <-annoBed[!((duplicated(annoBed$start)|duplicated(annoBed$start, fromLast=T)) & (duplicated(annoBed$end)|duplicated(annoBed$end, fromLast=T))),]
      annoBed$peakNumber <- c(paste0("Peak",1:length(annoBed$seqnames)))
      annoGR <- GRanges(seqnames=annoBed$seqnames, ranges=IRanges(start = annoBed$start, end = annoBed$end),
                        mcols=annoBed[grep("score", colnames(annoBed)):grep("GENENAME", colnames(annoBed))])
      annoGR <- annoGR[annoGR$mcols.TotalSignal >= quantile(annoGR$mcols.TotalSignal, probs=quant)]
      dists <- distanceToNearest(annoGR, EHGR)
      annoGR$mcols.distanceToTSS <- abs(annoGR$mcols.distanceToTSS)
      annoGR$EnhancerID <- c(EHGR$name[subjectHits(dists)])
      annoGR$EnhancerRanges <- c(paste0(seqnames(EHGR)[subjectHits(dists)], ":", start(EHGR)[subjectHits(dists)],
                                          "-" ,end(EHGR)[subjectHits(dists)]))
      annoGR$EnhancerScore <- c(score(EHGR)[subjectHits(dists)])
      annoGR$distToEnhancer <- c(as.data.frame(dists)$distance)
      annoGR$distLevel[annoGR$distToEnhancer == 0] <- EH_level_order[1]
      annoGR$distLevel[annoGR$distToEnhancer <= 50 & annoGR$distToEnhancer > 0] <- EH_level_order[2]
      annoGR$distLevel[annoGR$distToEnhancer <= 200 & annoGR$distToEnhancer > 50] <- EH_level_order[3]
      annoGR$distLevel[annoGR$distToEnhancer <= 500 & annoGR$distToEnhancer > 200] <- EH_level_order[4]
      annoGR$distLevel[annoGR$distToEnhancer <= 1000 & annoGR$distToEnhancer > 500] <- EH_level_order[5]
      annoGR$distLevel[annoGR$distToEnhancer <= 2500 & annoGR$distToEnhancer > 1000] <- EH_level_order[6]
      annoGR$distLevel[annoGR$distToEnhancer <= 5000 & annoGR$distToEnhancer > 2500] <- EH_level_order[7]
      annoGR$distLevel[annoGR$distToEnhancer <= 10000 & annoGR$distToEnhancer > 5000] <- EH_level_order[8]
      annoGR$distLevel[annoGR$distToEnhancer > 10000] <- EH_level_order[9]
      annoGR$distLevel <- factor(annoGR$distLevel, levels = EH_level_order)
      annoGR$graphAnnotation <- annoGR$mcols.annotation
      annoGR$graphAnnotation[grepl("Intron",annoGR$graphAnnotation)] <- "Intronic"
      annoGR$graphAnnotation[grepl("Exon", annoGR$graphAnnotation)] <- "Exonic"
      
      annoGR$distLevel <- factor(annoGR$distLevel, levels = EH_level_order)
      annoGR$graphAnnotation <- factor(annoGR$graphAnnotation, levels = level_order)
      print(paste0("Finished annotating ", histName, " ", quant, " ", enhancerName))
      {
      { # UCSC bar chart
        png(file=paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName, "/", histName, "_annotation_UCSC_knowGene_", quant, "quantile_barChart.png"), width = 1920,
            height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
        
        barChart <- ggplot() + geom_bar(aes(x=1, fill=factor(annoGR$graphAnnotation, level = rev(level_order))), position = "stack")+
          scale_y_continuous(breaks = cumsum(summary(factor(annoGR$graphAnnotation, level = level_order)))
                             [summary(factor(annoGR$graphAnnotation, level = level_order)) >
                                 quantile(summary(factor(annoGR$graphAnnotation, level = rev(level_order))), probs= 0.75)],
                             labels = round(cumsum(summary(factor(annoGR$graphAnnotation, level = level_order)))/length(annoGR$graphAnnotation), digits = 2)
                             [summary(factor(annoGR$graphAnnotation, level = level_order)) >
                                 quantile(summary(factor(annoGR$graphAnnotation, level = rev(level_order))), probs= 0.75)]) +
          scale_fill_manual(labels = rev(level_order), values=colourVec)+
          theme_bw()+
          theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                #axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                plot.title = element_text(hjust = 0.5)
          )+
          labs(title=paste0("Annotation Type Assignment For ", histName, " Peaks to UCSC KnownGene Database (Quantile = ", quant, ")"), y = "Proportion of Annotation Type",
               fill ="Annotation Type")+
          guides(fill=guide_legend(reverse = TRUE))+
          coord_flip()
        
        print(barChart + theme_bw(base_size = 20))
        dev.off()
        print(paste0("Printing UCSC BarChart ", histName, " ", quant, " ", enhancerName ))
      }
    
    { # EH BarChart
      png(file=paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName, "/", histName, "_annotation_enhancerAtlas_", enhancerName,"_", quant, "quantile_barChart.png"), width = 1920,
          height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
      
      EHbarChart <- ggplot() + geom_bar(aes(x=1, fill=factor(annoGR$distLevel, level = rev(EH_level_order))), position = "stack")+
        scale_y_continuous(breaks = cumsum(summary(factor(annoGR$distLevel, level = EH_level_order)))
                           [summary(factor(annoGR$distLevel, level = EH_level_order)) >
                               quantile(summary(factor(annoGR$distLevel, level = rev(EH_level_order))), probs= 0.75)],
                           labels = round(cumsum(summary(factor(annoGR$distLevel, level = EH_level_order)))/length(annoGR$distLevel), digits = 2)
                           [summary(factor(annoGR$distLevel, level = EH_level_order)) >
                               quantile(summary(factor(annoGR$distLevel, level = rev(EH_level_order))), probs= 0.75)]) +
        scale_fill_manual(labels = rev(EH_level_order), values=colourVec)+
        theme_bw()+
        theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
              #axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.title = element_text(hjust = 0.5)
        )+
        labs(title=paste0("Annotation Type Assignment For ", histName, " Peaks To enhancerAtlas Enhancers (", enhancerName, ", Quantile = ", quant, ")"), y = "Proportion of Annotation Type",
             fill ="Annotation Type")+
        guides(fill=guide_legend(reverse = TRUE))+
        coord_flip()
      
      print(EHbarChart + theme_bw(base_size = 20))
      dev.off()
      print(paste0("Printing EH BarChart ", histName, " ", quant, " ", enhancerName ))
    }
    {# UCSC pie chart
      {
        png(file=paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName, "/", histName, "_annotation_UCSC_knowGene_", quant, "quantile_PieChart.png"), width = 1920,
            height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
        
        prop = rev(summary(factor(annoGR$graphAnnotation, level = level_order)))
        labs = rev(summary(factor(annoGR$graphAnnotation, level = level_order)))
        ypos = cumsum(prop)- prop*0.5
        textDF <- data.frame(ypos, labs)
        PieChart <- ggplot()+
          geom_bar(aes(x="", fill=factor(annoGR$graphAnnotation, level = level_order)))+
          labs(title=paste0("Assignment of Annotation Types For CUT&Tag ", histName, " Peaks To UCSC KnownGene Database (Quantile = ", quant, ")"), y = "Proportion of Annotation Type",
               fill ="Annotation Type")+
          scale_fill_manual(
            #labels =level_order,
            values=colourVec)+
          theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                plot.title = element_text(size=15, hjust = 0.5))+
          scale_y_continuous(labels = scales::percent_format(scale = 1))+
          geom_text(data=subset(textDF, textDF$labs>100), aes(x =1.6,y=ypos, label=labs), color="black", size=8)+
          theme_void()+
          coord_polar("y", start=0)
        
        print(PieChart + theme_void(base_size = 15))
        dev.off()
        print(paste0("Printing UCSC Piechart ", histName, " ", quant, " ", enhancerName ))
      }
    }
    { # EH pie
      png(file=paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName, "/", histName, "_annotation_enhancerAtlas_", enhancerName,"_", quant, "quantile_PieChart.png"), width = 1920,
          height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
      
      prop = rev(summary(factor(annoGR$distLevel, level = EH_level_order)))
      labs = rev(summary(factor(annoGR$distLevel, level = EH_level_order)))
      ypos = cumsum(prop)- prop*0.5
      textDF <- data.frame(ypos, labs)
      EHPieChart <- ggplot()+
        geom_bar(aes(x="", fill=factor(annoGR$distLevel, level = EH_level_order)))+
        labs(title=paste0("Assignment of Annotation Types For CUT&Tag ", histName, " Peaks To enhancerAtlas Enhancers (", enhancerName, ", Quantile = ", quant, ")"), y = "Proportion of Annotation Type",
             fill ="Annotation Type")+
        scale_fill_manual(labels =EH_level_order, values=colourVec)+
        theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
              plot.title = element_text(size=15, hjust = 0.5))+
        scale_y_continuous(labels = scales::percent_format(scale = 1))+
        geom_text(data=subset(textDF, textDF$labs>100), aes(x =1.6,y=ypos, label=labs), color="black", size=8)+
        theme_void()+
        coord_polar("y", start=0)
      
      print(EHPieChart + theme_void(base_size = 15))
      dev.off()
      print(paste0("Printing EH Piechart ", histName, " ", quant, " ", enhancerName ))
    }
    { #sankey
      png(file=paste0(projPath, "/plots/Annotation/enhancerAtlas/", histName, "/", histName, "_annotation_UCSC_enhancerAtlas_", enhancerName,"_", quant, "quantile_sankeyChart.png"), width = 1920,
          height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
      
      DF<- data.frame("Sample" = c(rep(histName, length(annoGR[!(is.na(annoGR$EnhancerID))]))),
                      "UCSC_KnownGene_Annotation"=factor(annoGR$graphAnnotation, levels = level_order),
                      "Enhancer_Annotation"=factor(annoGR$distLevel, levels = EH_level_order))
      sankDF <- DF %>% make_long(Sample, UCSC_KnownGene_Annotation, Enhancer_Annotation)
      dagg <- sankDF%>% dplyr::group_by(node)%>%tally()
      sankDF2 <- merge(sankDF, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
      sankDF2$node <- factor(sankDF2$node, levels= rev(c(histName, level_order, EH_level_order)))
      sankDF2$next_node <- factor(sankDF2$next_node, levels= rev(c(histName, level_order, EH_level_order)))
      SankeyPlot <- ggplot(sankDF2, aes(x = x,
                                        next_x = next_x,
                                        node = node,
                                        next_node = next_node,
                                        fill = factor(node),
                                        label=paste0(node," n=", n))) +
        geom_sankey(flow.alpha = 0.5, node.color = 1) +
        labs(title=paste0("Assignment of Annotation Types For ", histName," CUT&Tag To UCSC KnownGene And enhancerAtlas (", enhancerName, ", Quantile = ", quant, ")"))+
        scale_fill_manual(values =colourVec,breaks=c(level_order, EH_level_order), guide=guide_legend(title = "Annotation", ncol=2, nrow=13, direction="vertical"))+
        geom_sankey_label(size = 5, color = 1, fill = "white") +
        theme_sankey(base_size=16)+
        theme(axis.title.x = element_blank(),
              plot.title = element_text(size=20))
      
      print(SankeyPlot + theme_void(base_size = 20))
      dev.off()
      print(paste0("Printing Sankey chart ", histName, " ", quant, " ", enhancerName ))
    }
  }
  
  # so we got the graphs
  # just get the interactions GRlist and map enhancers to enhancer interactions
  if(quant==0)
  {
    if(enhancerName != "IPSC")
    {
    preEHint <- read.table(file = enhancerInteractionFileNames[grepl("H1", enhancerInteractionFileNames)], sep="\t", quote="")
    preEHint$EnhancerRanges <- unlist(strsplit(preEHint$V1, "_"))[grepl(":",unlist(strsplit(preEHint$V1, "_")))]
    preEHint$interactedGENENAME <- unlist(strsplit(preEHint$V1, "\\$"))[c(FALSE,TRUE,FALSE,FALSE,FALSE)]
    preEHint$interactedTSSposition <- unlist(strsplit(preEHint$V1, "\\$"))[c(FALSE,FALSE,FALSE,TRUE,FALSE)]
    preEHint$interactionScore <- preEHint$V2
    EHintGR <- as(preEHint[,3:6]$EnhancerRanges, "GRanges")
    EHintGR$interactedGeneName <- preEHint[,3:6]$interactedGENENAME
    EHintGR$interactedTSSposition <- preEHint[,3:6]$interactedTSSposition
    EHintGR$interactionScore <- preEHint[,3:6]$interactionScore
    EHintGR <- subsetByOverlaps(EHintGR, blacklist.hg38, type="any", invert=TRUE)
    interactionHits <- findOverlaps(EHintGR, EHGR, type = "any")
    EHintGR$EnhancerID <- NA
    EHintGR[queryHits(interactionHits)]$EnhancerID <- EHGR$name[subjectHits(interactionHits)]
    
    testDF <-data.frame(interactedGeneName = EHintGR$EnhancerID,
                        interactionScore = EHintGR$interactionScore,
                        interactedTSSposition = EHintGR$interactedTSSposition,
                        interactedGeneName = EHintGR$interactedGeneName,
                        EnhancerID = EHintGR$EnhancerID)
    newlyAnnoGR <- merge(annoGR, testDF, by.x="EnhancerID", by.y="EnhancerID")
    newlyAnnoGR
    
    if(!dir.exists(paste0(dirname(fileName), "/enhancerAnnotated/EAaligned/base/")))
    {
      mkdirs(paste0(dirname(fileName), "/enhancerAnnotated/EAaligned/base/"))
    } 
    write.table(newlyAnnoGR, file= paste0(dirname(fileName), "/enhancerAnnotated/EAaligned/base/", histName, "_", enhancerName, "_", quant, "_EAaligned.txt"),sep="\t", row.names = F, col.names=T, quote=FALSE)
    }else{write.table(annoGR, file= paste0(dirname(fileName), "/enhancerAnnotated/EAaligned/base/", histName, "_", enhancerName, "_", quant, "_EAaligned.txt"),sep="\t", row.names = F, col.names=T, quote=FALSE)}
  }
  print(paste0("Quant was ", quant))
  }
  }
  print(paste0("HistName was ", histName))
}


  

