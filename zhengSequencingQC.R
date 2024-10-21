# R Code Adapted from Zheng C&T Tutorial From Protocols.io
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis) # for colours
library(GenomicRanges)
library(ggpubr) ## For customizing figures

#### Alignment Summary
projPath="C:/Users/mjmka/Documents/Bioinformatics"
sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX", "IPSC_2O", "D7_2O"
               #,"Undetermined"
               )
histList = c("H3K27me3", "H3K27ac", "yH2AX", "2O")

## Collect the alignment results from the bowtie2 alignment summary files
alignResult = c()
for(hist in sampleList){
  alignRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[2], Cell = histInfo[1], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))
write.csv(alignResult, file=paste0(projPath, "/plots/SequencingQA/AlignmentResult.csv"))

# Spike-In Alignment
# spikeAlign = c()
# for(hist in sampleList){
#   spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"), header = FALSE, fill = TRUE)
#   alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
#   histInfo = strsplit(hist, "_")[[1]]
#   spikeAlign = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
#                           SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
#                           MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
#                           AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
# }
# spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
# # spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
# ##=== R command ===## 
# alignSummary = left_join(alignResult, spikeAlign, by = c("Histone", "Replicate", "SequencingDepth")) %>%
#   mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"), 
#          AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
# alignSummary

## Generate sequencing depth boxplot
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 35) +
  ylab("Sequencing Depth per Million") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38/1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 35) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (hg38)")

fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 35) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (hg38)")

# fig3D = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
#   geom_boxplot() +
#   geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 18) +
#   ylab("Spike-in Alignment Rate") +
#   xlab("") +
#   ggtitle("D. Alignment Rate (E.coli)")

png(file=paste0(projPath, "/plots/SequencingQA/", "SequencingQA_NoUndetermined_SeqDepth.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")

print(fig3A)
dev.off()

png(file=paste0(projPath, "/plots/SequencingQA/", "SequencingQA_NoUndetermined_MappedFragNum.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")

print(fig3B)
dev.off()

png(file=paste0(projPath, "/plots/SequencingQA/", "SequencingQA_NoUndetermined_AlignmentRate.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")

print(fig3C)
dev.off()

print(ggarrange(fig3A, fig3B, fig3C, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom"))
## Summarize the duplication information from the picard summary outputs.
dupResult = c()
for(hist in sampleList){
  dupRes = read.table(paste0(projPath, "/alignment/rmDuplicate/picard_summary/", hist, "_picard.rmDup.txt"), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[2], Cell = histInfo[1], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% 
                           as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% 
                           as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>%
    mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}
dupResult$Histone = factor(dupResult$Histone, levels = histList)
alignDupSummary = left_join(alignResult, dupResult, by = c("Histone", "Cell", "MappedFragNum_hg38")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))
alignDupSummary
write.csv(alignDupSummary, file=paste0(projPath, "/plots/SequencingQA/DuplicationResult.csv"))



## generate boxplot figure for the  duplication rate
fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 35) +
  ylab("Duplication Rate (*100%)") +
  xlab("") 

fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 35) +
  ylab("Estimated Library Size") +
  xlab("") 

fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 35) +
  ylab("# of Unique Fragments") +
  xlab("")

png(file=paste0(projPath, "/plots/SequencingQA/", "SequencingQA_NoUndetermined_DuplicationRate.png"), width = 960,
    height=540,  pointsize = 35,antialias = "subpixel", type = "cairo-png")

ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
dev.off()



## Collect the fragment size information
fragLen = c()
for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.rmDup.txt"), header = FALSE) %>%
    mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[2], Cell = histInfo[1], sampleInfo = hist) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)
## Generate the fragment size density plot (violin plot)
fig5A = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Cell)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

png(file=paste0(projPath, "/plots/SequencingQA/", "SequencingQA_NoUndetermined_FragmentLength.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
ggarrange(fig5A, fig5B, ncol = 2)
dev.off()

#### we don't have any replicates so we don't need this
# library(corrplot) ## For correlation plot
# for reproducibility between replicates
# reprod = c()
# fragCount = NULL
# for(hist in sampleList){
# 
#   if(is.null(fragCount)){
# 
#     fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
#     colnames(fragCount) = c("chrom", "bin", hist)
# 
#   }else{
# 
#     fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
#     colnames(fragCountTmp) = c("chrom", "bin", hist)
#     fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
# 
#   }
# }
# 
# M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 
# 
# corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

# ##=== R command ===## 
# scaleFactor = c()
# multiplier = 10000
# for(hist in sampleList){
#   spikeDepth = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.seqDepth"), header = FALSE, fill = TRUE)$V1[1]
#   
#   histInfo = strsplit(hist, "_")[[1]]
#   scaleFactor = data.frame(scaleFactor = multiplier/spikeDepth, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(scaleFactor, .)
# }
# scaleFactor$Histone = factor(scaleFactor$Histone, levels = histList)
# left_join(alignDupSummary, scaleFactor, by = c("Histone", "Replicate"))

##=== R command ===##
## Generate sequencing depth boxplot
# fig6A = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
#   geom_boxplot() +
#   geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 20) +
#   ylab("Spike-in Scalling Factor") +
#   xlab("")
# 
# normDepth = inner_join(scaleFactor, alignResult, by = c("Histone", "Replicate")) %>% mutate(normDepth = MappedFragNum_hg38 * scaleFactor)
# 
# fig6B = normDepth %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
#   geom_boxplot() +
#   geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 20) +
#   ylab("Normalization Fragment Count") +
#   xlab("") + 
#   coord_cartesian(ylim = c(1000000, 130000000))
# ggarrange(fig6A, fig6B, ncol = 2, common.legend = TRUE, legend="bottom")

##=== R command ===## 
peakN = c()
peakWidth = c()
peakType = c(
 # "control", # we don't use the control seacr peaks
  "top0.01")
peakNature=c("relaxed", "stringent")
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  if(histInfo[2] != "2O"){
    for(type in peakType){
      for(nature in peakNature){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/CUTandTag/", hist, "_seacr_", type, "_", nature,".peaks.", nature,".bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakNature = nature, peakType = type, Histone = histInfo[2], Cell = histInfo[1]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakNature = nature ,peakType = type, Histone = histInfo[2], Cell = histInfo[1])  %>% rbind(peakWidth, .)
      }
    }
  }
}
peakN %>% select(Histone, Cell, peakNature, peakN)
#### we don't have replicates so we don't use this
##=== R command ===## 
# library(chromVAR)
# 
# bamDir = paste0(projPath, "/alignment/bam")
# inPeakData = c()
# ## overlap with bam file to get count
# for(hist in histL){
#   for(nature in peakNature){
#     peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
#     peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
#     bamFile = paste0(bamDir, "/", hist, "_", rep, "_bowtie2.mapped.bam")
#     fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
#     inPeakN = counts(fragment_counts)[,1] %>% sum
#     inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = hist, Replicate = rep))
#   }
#}

frip = left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% mutate(frip = inPeakN/MappedFragNum_hg38 * 100)
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)

fig7A = peakN %>% ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Cell), position = position_jitter(0.15)) +
  facet_grid(~peakNature) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

fig7B = peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Cell~peakNature) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

# fig7C = peakReprod %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
#   geom_bar(stat = "identity") +
#   geom_text(vjust = 0.1) +
#   facet_grid(Replicate~peakType) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 18) +
#   ylab("% of Peaks Reproduced") +
#   xlab("")

# fig7D = frip %>% ggplot(aes(x = Histone, y = frip, fill = Histone, label = round(frip, 2))) +
#   geom_boxplot() +
#   geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
#   scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
#   scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
#   theme_bw(base_size = 18) +
#   ylab("% of Fragments in Peaks") +
#   xlab("")

png(file=paste0(projPath, "/plots/SequencingQA/", "SequencingQA_NoUndetermined_PeakStats.png"), width = 1920,
    height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
ggarrange(fig7A, fig7B, 
         # fig7C, fig7D, ncol = 2, nrow=2,
         ncol=2, nrow=1,
          common.legend = TRUE, legend="bottom")
while (!is.null(dev.list()))  dev.off()

### just some code I used to plot score distributions
# source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))
# fileNames <- findInDirList(paste0(projPath, "/peakCalling/SEACR"), "CUTandTag")
# fileNames <- fileNames[!(grepl("stringent|2O|Undetermined", fileNames))]
# fileNames
# for (fileName in fileNames)
# {
#   #read data
#   histName <- sampleList[str_detect(fileName, sampleList)]
#   DF <- read.table(file = fileName, header = F, sep="\t", quote="")
#   scoreDF <- data.frame(score = sort(DF$V4), SampleName = c(histName))
#   scoreDensity <- ggplot(scoreDF, mapping = aes(x= score))+
#     geom_density()+
#     scale_x_log10()
#   png(file=paste0(projPath, "/plots/SequencingQA/",histName,"_score_distribution.png"), width = 1920,
#       height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
#   print(scoreDensity)
#   dev.off()
# }
