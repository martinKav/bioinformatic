# SEdb Preparation and Investigation
# Author: Martin Kavanagh
# Date: 06/08/2024
# ******************************************
library(stringr)
library(readr)
library(ggplot2)
library(scales)
library(viridis)
library(ggblend)
library(ggdist)
library(wesanderson)
# read SEdb files
projPath="C:/Users/mjmka/Documents/Bioinformatics"
sampleList=c("IPSC_H3K27ac", "IPSC_H3K27me3", "IPSC_yH2AX"
             #, "IPSC_2O"
             , "D7_H3K27ac", "D7_H3K27me3", "D7_yH2AX"
             #, "D7_2O", "Undetermined"
             )
EnhancerTypeList=c("Enhancer"
                   #, "PromoterEnhancer"
                   )
plusminusList=(c("15000", "10000", "5000", "2000", "1000", "500"
                 #,"0"
                 ))
for (enhancerType in EnhancerTypeList)
{
  for (plusminus in plusminusList)
  {
     # load data
    histDF <- data.frame()
    for (histName in sampleList)
    {
      tempDF <- read.table(file=paste0(projPath, "/data/Genehancer/plots/", histName,"_genehancer_",enhancerType,"_plusminus", plusminus,"bp_start_extended.window.counts.txt"))
      tempDF$name <- histName
      histDF <- rbind(histDF, tempDF)
    }
    histDF$name <- as.factor(histDF$name)
    
    
    png(file=paste0(projPath, "/plots/Genehancer/", "Genehancer_", enhancerType, "_", plusminus,"bp_HistDistribution.png"), width = 1920,
        height=1080,  pointsize = 35,antialias = "subpixel", type = "cairo-png")
    plot <- ggplot(histDF, aes(x = V1, y= V2, color = name))+
      ggtitle(paste0("Plot of Histone Distribution at ", enhancerType,"s from Genehancer Dataset")) +
      labs(x ="Distribution of Histone at Enhancer", y = "Count", color = "Sample")+
      scale_color_manual(values = c("red", "orange", "blue", "green", "pink", "purple", "cyan", "magenta", "brown"))+
      geom_point(size = 2, alpha = 0.9) * (blend("lighten") + blend("multiply", alpha = 0.5))+ 
      guides(colour = guide_legend(override.aes = list(size=10)))+
      theme_gray(base_size=12*(1920/1080))
      #geom_smooth()
    print(plot,)
    dev.off()
  }
}

avg_enhancerPlot <- ggplot()+
  ggtitle("Plot of Histone Distribution at Enhancers from Genehancer Dataset") +
  labs(x ="Distribution of Histone at Enhancer", y = "Count", colour = "sample")+
  geom_point(size = 0.1, alpha = 0.1)
