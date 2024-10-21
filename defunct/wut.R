# R script for enriching annotated gene sets
# enrich genesets

library(dplyr)
library(stringr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(rentrez)
library(rGREAT)
library(biomaRt)
library(data.table)
library(ggplot2)
projPath="C:/Users/mjmka/Documents/Bioinformatics"
source(paste0(projPath, "/scripts/subscripts/findInDirList.R"))

# can do great, kegg, meme
# clustering

# we cant infer on gene expression from this I dont think - because the ensemblIDs wont match the symbols gotten from the interactions
# a main limitations is that we can only take the first ensemblID given to avoid duplication - may not be accurate but difficult to determine via cis regulation
# find data
fileNames <- findInDirList(paste0(projPath, "/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/base"), ".txt")


for (fileName in fileNames)
{
  #read data
  DF <- read.table(file = fileName, header = T, sep="\t", quote="")
  #insert to track duplication
  DF[duplicated(DF$mcols.peakNumber)|duplicated(DF$mcols.peakNumber, fromLast=T),]
  geneIDs2 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= DF$interactedGeneName, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))  
  geneIDs2 <- geneIDs2[!(grepl("LRG", geneIDs2$GENEID)),]
  geneIDs2 <-geneIDs2[!(duplicated(geneIDs2$SYMBOL)),]
  # need to de-duplicate the ensembl IDs or else they will bias the downstream gene ontology
  colnames(geneIDs2) <- c("SYMBOL", "ensemblID")
  entrezIDs = mapIds(org.Hs.eg.db,
                       keys=geneIDs2$ensemblID, #Column containing Ensembl gene ids
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
  entrezDF <- data.frame( ensemblID=names(entrezIDs),entrezID =entrezIDs, row.names = NULL)
  addedIDs <- merge(geneIDs2, entrezDF, by="ensemblID" )
  colnames(addedIDs) <- c("interactedGeneENSEMBL", "interactedGeneName", "interactedGeneENTREZ")
  DF2 <-merge(DF, addedIDs, by="interactedGeneName")
  write.table(DF2,file=paste0(projPath, "/annotationsAndEnrichment/SEACR/CUTandTag/enhancerAnnotated/GHaligned/EnsemblEntrezAnnotated/",strsplit(basename(fileName), ".txt")[[1]], "_EnsemblEntrez.txt"),sep="\t", row.names = F, col.names=T, quote=FALSE)
}

