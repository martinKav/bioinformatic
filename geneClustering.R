# biocor cluster analysis
library(dplyr)
library(stringr)
library(readr)
library(GenomicRanges)
library(BioCor)

#prepare sample - try these first
c(
  "~/Bioinformatics/geneOverlaps/IPSC_H3K27ac-D7_H3K27ac-IPSC_yH2AX/IPSC_H3K27acInIPSC_yH2AX_withoutD7_H3K27ac_overlap_quant0.99.txt",
  "~/Bioinformatics/geneOverlaps/IPSC_H3K27ac-D7_H3K27ac-IPSC_yH2AX/IPSC_H3K27acInIPSC_yH2AX_withoutD7_H3K27ac_overlap_quant0.txt",
  "~/Bioinformatics/geneOverlaps/IPSC_H3K27ac-D7_H3K27ac-IPSC_yH2AX/IPSC_H3K27acInIPSC_yH2AX_withoutD7_H3K27ac_overlap_quant0.5.txt",
  "~/Bioinformatics/geneOverlaps/IPSC_H3K27ac-D7_H3K27ac-IPSC_yH2AX/IPSC_H3K27acInIPSC_yH2AX_withoutD7_H3K27ac_overlap_quant0.9.txt",
  "~/Bioinformatics/geneOverlaps/IPSC_H3K27ac-D7_H3K27ac-IPSC_yH2AX/IPSC_H3K27acInIPSC_yH2AX_withoutD7_H3K27ac_overlap_quant0.75.txt"
)
entrezIDs <- read.table(file="~/Bioinformatics/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/
                        GHaligned/EnsemblEntrezAnnotated/IPSC_yH2AX_EnhancersAndPromoterEnhancers_0_GHaligned_EnsemblEntrez.txt", sep="\t", header=T, quote="" )
bed <- read.table(file="~/Bioinformatics/geneOverlaps/IPSC_H3K27ac-D7_H3K27ac-IPSC_yH2AX/IPSC_H3K27acInIPSC_yH2AX_withoutD7_H3K27ac_overlap_quant0.99.txt", sep="\t", header=T, quote="" )
entrezIDs[entrezIDs$interactedGeneName %in% bed$GENENAMES]

tempDF <- data.frame(entrezIDs = entrezIDs$interactedGeneENTREZ[unique(entrezIDs$interactedGeneENTREZ[entrezIDs$interactedGeneName %in% bed$GENENAMES])],
                     GENENAMES <- entrezIDs[unique(entrezIDs$interactedGeneENTREZ[entrezIDs$interactedGeneName %in% bed$GENENAMES]),]$interactedGeneName)
tempDF[!is.na(tempDF$entrezIDs),]
newBed <- merge(tempDF, bed)
newBed
unique(newBed$entrezIDs[!is.na(newBed$entrezIDs)])
entrezIDs$interactedGeneENTREZ[unique(entrezIDs$interactedGeneENTREZ[entrezIDs$interactedGeneName %in% bed$GENENAMES])]
                     entrezIDs[unique(entrezIDs$interactedGeneENTREZ[entrezIDs$interactedGeneName %in% bed$GENENAMES]),]$interactedGeneName


entrezList <- unique(newBed$entrezIDs[!is.na(newBed$entrezIDs)])





## Load libraries with the data of the pathways
library("org.Hs.eg.db")
library("reactome.db")
genesKegg <- as.list(org.Hs.egPATH)
genesReact <- as.list(reactomeEXTID2PATHID)
genesReact

genesReact <- lapply(genesReact, function(x){
  unique(grep("R-HSA-", x, value = TRUE))
})
genesReact <- genesReact[lengths(genesReact) >= 1] 


geneSim("672", "675", genesKegg)
## [1] 0.0824295
geneSim("672", "675", genesReact)

mgeneSim(c("BRCA1" = "672", "BRCA2" = "675", "NAT2" = "10"), genesKegg)

mgeneSim(c("BRCA1" = "672", "BRCA2" = "675", "NAT2" = "10"), genesReact)


# clustering example - ok this is pretty portable
genes.id <- as.character(entrezList)
genes.id <- mapIds(org.Hs.eg.db, keys = genes.id, keytype = "ENTREZID", 
                   column = "SYMBOL")
## 'select()' returned 1:1 mapping between keys and columns
genes <- names(genes.id)
names(genes) <- genes.id
react <- mgeneSim(genes, genesReact)
kegg <- mgeneSim(genes, genesKegg)
## We remove genes which are not in list (hence the warning):
nan <- genes %in% names(genesReact)
react <- react[nan, nan]
hc <- hclust(as.dist(1 - react))

plot(hc, main = "Similarities between genes - REACT") # shows relationship between genes
# cluster anal
mycl <- cutree(hc, h = 0.2)
mycl
clusters <- split(genes[nan], as.factor(mycl)) # puts them into a list
clusters
(clusters <- clusters[lengths(clusters) >= 2]) # removes clusters with only 1 gene
names(clusters) <- paste0("cluster", names(clusters)) #names them cluster1:8
clusters
sim_clus1 <- mclusterSim(clusters, genesReact) # then clusters the clusters
plot(hclust(as.dist(1 - sim_clus1)), 
     main = "Similarities between clusters by pathways") # shows similarities between clusters by pathways
sim_clus2 <- mclusterGeneSim(clusters, genesReact)
plot(hclust(as.dist(1 - sim_clus2)), 
     main ="Similarities between clusters by genes") #shows similarities between clusters by genes

kegg <- kegg[rowSums(is.na(kegg)) != ncol(kegg), rowSums(is.na(kegg)) != ncol(kegg)]
colnames(kegg)
hc <- hclust(as.dist(1-kegg))
plot(hc, main = "Similarlities between genes - KEGG")
mycl <- cutree(hc, h = 0.2)
mycl
clusters <- split(genes[colnames(kegg)], as.factor(mycl)) # puts them into a list
(clusters <- clusters[lengths(clusters) >= 2]) # removes clusters with only 1 gene
names(clusters) <- paste0("cluster", names(clusters)) #names them cluster1:8
sim_clus1 <- mclusterSim(clusters, genesKegg) # then clusters the clusters
plot(hclust(as.dist(1 - sim_clus1)), 
     main = "Similarities between clusters by pathways") # shows similarities between clusters by pathways
sim_clus2 <- mclusterGeneSim(clusters, genesKegg)
plot(hclust(as.dist(1 - sim_clus2)), 
     main ="Similarities between clusters by genes") #shows similarities between clusters by genes
