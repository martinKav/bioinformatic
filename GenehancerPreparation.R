# Genehancer & SEdb & EnhancerAtlas Preparation
# Author: Martin Kavanagh
# Date: 06/08/2024
# **************************************************************
#loading libs
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
library(GenomicRanges)
library(ggsankey)
library(R.utils)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
colourVec=c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "#FB9A99", "palegreen2", "cyan",  "#FDBF6F",  "peachpuff3", 
            "khaki2", "maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")
GHtypeList=c("EnhancerOnly", "EnhancersAndPromoterEnhancers")
quantList=c(0, 0.5, 0.75, 0.9, 0.99)
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX",
               "H1_H3K27ac","H1_H3K27me3", "H1_H3K4me1", "H9_H3K27ac","H9_H3K27me3", "H9_H3K4me1", "IPSC_ORLANDO_yH2AX")
level_order <- c('Promoter', 'Promoter (<=1kb)', 'Promoter (1-2kb)', 'Promoter (2-3kb)', "Promoter (3-4kb)", "Promoter (4-5kb)",
                 "Promoter (5-6kb)", "Promoter (6-7kb)", "3' UTR", "Exonic", "Intronic", "5' UTR", "Downstream (<=300bp)",
                 "Distal Intergenic") 
GH_level_order <- c("At GH Enhancer", ">50 bp from GH Enhancer","50-100 bp from GH Enhancer","200-500 bp from GH Enhancer","500-1k bp from GH Enhancer",
                    "1k-2.5k bp from GH Enhancer","2.5k-5k bp from GH Enhancer", "5k-10k bp from GH Enhancer", "<10k bp from GH Enhancer" )

# loading blacklist to blacklist genehancers inside the blacklisted regions
blackBed <- read.table(file = paste0(projPath, "/data/ncbi/hg38-blacklist.v2.bed" ), fill=T, header = F, sep = "\t")
colnames(blackBed) <- c("chrom", "start", "end", "reason")
blacklist.hg38 <- GRanges(ranges=IRanges(start=blackBed$start, end=blackBed$end),
                          seqnames=blackBed$chrom,
                          reason=blackBed$reason)


# loading the genehancer-gene interactions
GHint <- read.table( paste0(projPath, "/data/Genehancer/genehancer_interactions.txt"), header = TRUE, sep = "\t")
GHintGR <- GRanges(seqnames=GHint$chrom, IRanges(start=GHint$chromStart, end=GHint$chromEnd),
                   mcols=GHint[grep("name", colnames(GHint)):grep("geneStrand", colnames(GHint))])
GHintGR <- subsetByOverlaps(GHintGR, blacklist.hg38, type="any", invert=TRUE)

fileNames <- findInDirList(paste0(projPath, "/annotationsAndEnrichment/"), c("encode", "orlando", "CUTandTag"))
fileNames <- fileNames[!(grepl("_2O|Undetermined", fileNames))]
fileNames <- fileNames[!(grepl("enrichment", fileNames))]
fileNames <- fileNames[!(grepl("stringent", fileNames))]
fileNames <- fileNames[!(grepl("enhancerAnnotated", fileNames))]
fileNames
for (GHtype in GHtypeList)
{
  GHbed <- read.table(paste0(projPath, "/data/Genehancer/genehancer.bed"))
  GHbed <- GHbed[GHbed$V11 != "Promoter",]
  if ( GHtype == "EnhancerOnly"){GHbed <- GHbed[GHbed$V11 != "Promoter/Enhancer",]}
  GHGR <- GRanges(seqnames=GHbed$V1 ,ranges = IRanges(start=GHbed$V2, end= GHbed$V3),mcols= GHbed[,4:12])
  GHGR <- subsetByOverlaps(GHGR, blacklist.hg38, type = "any", invert = TRUE)
  for(fileName in fileNames)
  {
    for (quant in quantList)
    {
      histName <- sampleList[str_detect(fileName, sampleList)]
      if(!dir.exists(paste0(projPath, "/plots/Annotation/Genehancer/", histName)))
      {
        mkdirs(paste0(projPath, "/plots/Annotation/Genehancer/", histName))
      } 
      annoBed <- read.table(file=fileName, sep="\t", quote="")
      annoBed <-annoBed[!((duplicated(annoBed$start)|duplicated(annoBed$start, fromLast=T)) & (duplicated(annoBed$end)|duplicated(annoBed$end, fromLast=T))),]
      annoBed$peakNumber <- c(paste0("Peak",1:length(annoBed$seqnames)))
      annoGR <- GRanges(seqnames=annoBed$seqnames, ranges=IRanges(start = annoBed$start, end = annoBed$end),
                        mcols=annoBed[grep("score", colnames(annoBed)):grep("GENENAME", colnames(annoBed))])
      annoGR <- annoGR[seqnames(annoGR) != "chrY"]
    
      annoGR <- annoGR[annoGR$mcols.score >= quantile(annoGR$mcols.score, probs=quant)]
      dists <- distanceToNearest(annoGR, GHGR)
      annoGR$mcols.distanceToTSS <- abs(annoGR$mcols.distanceToTSS)
      annoGR$GenehancerID <- c(GHGR$mcols.V4[subjectHits(dists)], rep(NA, times =length(annoGR[seqnames(annoGR) == "chrY"])))
      annoGR$GenehancerRanges <- c(paste0(seqnames(GHGR)[subjectHits(dists)], ":", start(GHGR)[subjectHits(dists)],
                 "-" ,end(GHGR)[subjectHits(dists)]), rep(NA, times =length(annoGR[seqnames(annoGR) == "chrY"])))
      annoGR$GenehancerScore <- c(GHGR$mcols.V5[subjectHits(dists)], rep(NA, times =length(annoGR[seqnames(annoGR) == "chrY"])))
      annoGR$distToGenehancer <- c(as.data.frame(dists)$distance, rep(NA, times =length(annoGR[seqnames(annoGR) == "chrY"])))
      annoGR$GenehancerType <- c(GHGR$mcols.V11[subjectHits(dists)], rep(NA, times =length(annoGR[seqnames(annoGR) == "chrY"])))
      annoGR$distLevel[annoGR$distToGenehancer == 0] <- "At GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 50 & annoGR$distToGenehancer > 0] <- ">50 bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 200 & annoGR$distToGenehancer > 50] <- "50-100 bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 500 & annoGR$distToGenehancer > 200] <- "200-500 bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 1000 & annoGR$distToGenehancer > 500] <- "500-1k bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 2500 & annoGR$distToGenehancer > 1000] <- "1k-2.5k bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 5000 & annoGR$distToGenehancer > 2500] <- "2.5k-5k bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer <= 10000 & annoGR$distToGenehancer > 5000] <- "5k-10k bp from GH Enhancer"
      annoGR$distLevel[annoGR$distToGenehancer > 10000] <- "<10k bp from GH Enhancer"
      annoGR$distLevel <- factor(annoGR$distLevel, levels = GH_level_order)
      annoGR$graphAnnotation <- annoGR$mcols.annotation
      annoGR$graphAnnotation[grepl("Intron",annoGR$graphAnnotation)] <- "Intronic"
      annoGR$graphAnnotation[grepl("Exon", annoGR$graphAnnotation)] <- "Exonic"
      
      annoGR$distLevel <- factor(annoGR$distLevel, levels = GH_level_order)
      annoGR$graphAnnotation <- factor(annoGR$graphAnnotation, levels = level_order)
      print(paste0("Finished annotating ", histName, " ", quant, " ", GHtype))
      {
        if(GHtype=="EnhancerOnly")
        {
        { # UCSC bar chart
        png(file=paste0(projPath, "/plots/Annotation/Genehancer/", histName, "/", histName, "_annotation_UCSC_knowGene_", quant, "quantile_barChart.png"), width = 1920,
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
        print(paste0("Printing UCSC BarChart ", histName, " ", quant, " ", GHtype ))
        }
        }
        { # GH BarChart
        png(file=paste0(projPath, "/plots/Annotation/Genehancer/", histName, "/", histName, "_annotation_Genehancer_", GHtype,"_", quant, "quantile_barChart.png"), width = 1920,
            height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
  
        GHbarChart <- ggplot() + geom_bar(aes(x=1, fill=factor(annoGR$distLevel, level = rev(GH_level_order))), position = "stack")+
          scale_y_continuous(breaks = cumsum(summary(factor(annoGR$distLevel, level = GH_level_order)))
                             [summary(factor(annoGR$distLevel, level = GH_level_order)) >
                                 quantile(summary(factor(annoGR$distLevel, level = rev(GH_level_order))), probs= 0.75)],
                             labels = round(cumsum(summary(factor(annoGR$distLevel, level = GH_level_order)))/length(annoGR$distLevel), digits = 2)
                             [summary(factor(annoGR$distLevel, level = GH_level_order)) >
                                 quantile(summary(factor(annoGR$distLevel, level = rev(GH_level_order))), probs= 0.75)]) +
          scale_fill_manual(labels = rev(GH_level_order), values=colourVec)+
          theme_bw()+
          theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                #axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                plot.title = element_text(hjust = 0.5)
                )+
          labs(title=paste0("Annotation Type Assignment For ", histName, " Peaks To Genehancer Enhancers (", GHtype, ", Quantile = ", quant, ")"), y = "Proportion of Annotation Type",
               fill ="Annotation Type")+
          guides(fill=guide_legend(reverse = TRUE))+
          coord_flip()
  
        print(GHbarChart + theme_bw(base_size = 20))
        dev.off()
        print(paste0("Printing GH BarChart ", histName, " ", quant, " ", GHtype ))
        }
        {# UCSC pie chart
        if(GHtype == "EnhancerOnly")
        {
        png(file=paste0(projPath, "/plots/Annotation/Genehancer/", histName, "/", histName, "_annotation_UCSC_knowGene_", quant, "quantile_PieChart.png"), width = 1920,
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
        print(paste0("Printing UCSC Piechart ", histName, " ", quant, " ", GHtype ))
        }
        }
        { # GH pie
        png(file=paste0(projPath, "/plots/Annotation/Genehancer/", histName, "/", histName, "_annotation_Genehancer_", GHtype,"_", quant, "quantile_PieChart.png"), width = 1920,
            height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
  
        prop = rev(summary(factor(annoGR$distLevel, level = GH_level_order)))
        labs = rev(summary(factor(annoGR$distLevel, level = GH_level_order)))
        ypos = cumsum(prop)- prop*0.5
        textDF <- data.frame(ypos, labs)
        GHPieChart <- ggplot()+
          geom_bar(aes(x="", fill=factor(annoGR$distLevel, level = GH_level_order)))+
          labs(title=paste0("Assignment of Annotation Types For CUT&Tag ", histName, " Peaks To Genehancer Enhancers (", GHtype, ", Quantile = ", quant, ")"), y = "Proportion of Annotation Type",
               fill ="Annotation Type")+
          scale_fill_manual(labels =GH_level_order, values=colourVec)+
          theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),
                plot.title = element_text(size=15, hjust = 0.5))+
          scale_y_continuous(labels = scales::percent_format(scale = 1))+
          geom_text(data=subset(textDF, textDF$labs>100), aes(x =1.6,y=ypos, label=labs), color="black", size=8)+
          theme_void()+
          coord_polar("y", start=0)
  
        print(GHPieChart + theme_void(base_size = 15))
        dev.off()
        print(paste0("Printing GH Piechart ", histName, " ", quant, " ", GHtype ))
        }
        { #sankey
        png(file=paste0(projPath, "/plots/Annotation/Genehancer/", histName, "/", histName, "_annotation_UCSC_Genehancer_", GHtype,"_", quant, "quantile_sankeyChart.png"), width = 1920,
            height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
  
        DF<- data.frame("Sample" = c(rep(histName, length(annoGR[!(is.na(annoGR$GenehancerID))]))),
                        "UCSC_KnownGene_Annotation"=factor(annoGR$graphAnnotation, levels = level_order),
                        "Genehancer_Annotation"=factor(annoGR$distLevel, levels = GH_level_order))
        sankDF <- DF %>% make_long(Sample, UCSC_KnownGene_Annotation, Genehancer_Annotation)
        dagg <- sankDF%>% dplyr::group_by(node)%>%tally()
        sankDF2 <- merge(sankDF, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
        sankDF2$node <- factor(sankDF2$node, levels= rev(c(histName, level_order, GH_level_order)))
        sankDF2$next_node <- factor(sankDF2$next_node, levels= rev(c(histName, level_order, GH_level_order)))
        SankeyPlot <- ggplot(sankDF2, aes(x = x,
                       next_x = next_x,
                       node = node,
                       next_node = next_node,
                       fill = factor(node),
                       label=paste0(node," n=", n))) +
          geom_sankey(flow.alpha = 0.5, node.color = 1) +
          labs(title=paste0("Assignment of Annotation Types For ", histName," CUT&Tag To UCSC KnownGene And Genehancer (", GHtype, ", Quantile = ", quant, ")"))+
          scale_fill_manual(values =colourVec,breaks=c(level_order, GH_level_order), guide=guide_legend(title = "Annotation", ncol=2, nrow=13, direction="vertical"))+
          geom_sankey_label(size = 5, color = 1, fill = "white") +
          theme_sankey(base_size=16)+
          theme(axis.title.x = element_blank(),
                plot.title = element_text(size=20))
  
        print(SankeyPlot + theme_void(base_size = 20))
        dev.off()
        print(paste0("Printing Sankey chart ", histName, " ", quant, " ", GHtype ))
        }
      }
      
      # so we got the graphs
      # just get the GHintGR and map enhancers to enhancer interactions
      if(quant==0)
      {
      enhSubset <- annoGR[annoGR$distToGenehancer <= 200]
      
      GHintGR # want name, score, value? Identifier, geneName
      
      testDF <-data.frame(interactionName = GHintGR$mcols.name,
                          interactionScore = GHintGR$mcols.score,
                          interactionValue = GHintGR$mcols.value,
                          interactedGeneName = GHintGR$mcols.geneName,
                          genehancerID = GHintGR$mcols.geneHancerIdentifier)
      newlyAnnoGR <- merge(annoGR, testDF, by.x="GenehancerID", by.y="genehancerID")
      if(!dir.exists(paste0(dirname(fileName), "/enhancerAnnotated/GHaligned/base/")))
      {
        mkdirs(paste0(dirname(fileName), "/enhancerAnnotated/GHaligned/base/"))
      } 
      write.table(newlyAnnoGR, file= paste0(dirname(fileName), "/enhancerAnnotated/GHaligned/base/", histName, "_", GHtype, "_", quant, "_GHaligned.txt"),sep="\t", row.names = F, col.names=T, quote=FALSE)
        }
    }
    print(paste0("Quant was ", quant))
  }
  print(paste0("HistName was ", histName))
}




GHbed <- read.table(paste0(projPath, "/data/Genehancer/genehancer.bed"))
# find what the enhancer type dist is
summary(factor(GHbed$V11))
plotDist <- data.frame(ID = GHbed$V4,type = as.factor(GHbed$V11), count = )
plotDist$count[plotDist$type == "Enhancer"] <- summary(factor(GHbed$V11))[1]
plotDist$count[plotDist$type == "Promoter"] <- summary(factor(GHbed$V11))[2]
plotDist$count[plotDist$type == "Promoter/Enhancer"] <- summary(factor(GHbed$V11))[3]
GHplot <- ggplot(plotDist, aes(x=type, y = count, fill=type))+
  labs(x="Genehancer Enhancer Classification", y = "Number of Enhancers of a Given Type",
       title = "Distribution of Different Genehancer Types")+
  geom_bar(position = 'dodge', stat='identity') +
  geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.25)+
  scale_fill_manual(values =c("red", "brown", "orange"))+
  theme_bw()

png(file=paste0(projPath, "/plots/Annotation/GenehancerDistributionBarChart.png"), width = 960,
    height=540,  pointsize = 20,antialias = "subpixel", type = "cairo-png")
print(GHplot + theme_bw(base_size = 20))
dev.off()
