# this script will take a GRlist, and return a GRlist with all the possible intersections of sets
library(stringr)
library(readr)
library(plyr)
library(dplyr)
library(GenomicRanges)
library(gtools)
source(paste0(projPath, "/scripts/subscripts/smallTools.R"))
getGRListIntersections <- function(newGRlist, inQuant=0, outQuant=0)
{
  # but how do you iterate the in and out lists
  holdList <- list()
  InOutList <- list()
  for (N in 1:length(names(newGRlist)))
  {
    x = names(newGRlist)
    x1 = combn(x, N)
    x2 <- matrix(ncol=NCOL(x1), nrow=(length(x) - NROW(x1)))
    if(length(x2) != 0)
    {
      for (col in 1:NCOL(x1)){ x2[,col] <- x[!(x %in% x1[,col])]}
    }
    list <- list(InList = x1, OutList = x2)
    InOutList[[N]] <- list
  }
  for (N in 1:length(names(newGRlist)))
  {
    for(col in 1:NCOL(InOutList[[N]][["InList"]]))
    {
      if(NCOL(InOutList[[N]][["InList"]]) > 1){InList <- InOutList[[N]][["InList"]][,col]}else{InList <- InOutList[[N]][["InList"]]}
      if(NCOL(InOutList[[N]][["OutList"]]) > 1){OutList <- InOutList[[N]][["OutList"]][,col]}else{OutList <- InOutList[[N]][["OutList"]]}
      intersectGR <- filterByTotalSignal(newGRlist[[InList[1]]], inQuant)
      
      for(inGRName in InList)
      {
        intersectGR <- subsetByOverlaps(intersectGR, filterByTotalSignal(newGRlist[[inGRName]],inQuant))
        
        # from the subject GR, we want to grab the total signal for each overlap and display it
        # but this was deemed too difficult, so it will just be easier to add it in later
        
      }
      if(length(OutList))
      {
        for(outGRName in OutList)
        {
          intersectGR <- subsetByOverlaps(intersectGR, filterByTotalSignal(newGRlist[[outGRName]], outQuant), invert = TRUE)
        }
      }
      
      if(length(OutList)){holdList[[paste0(paste0(InList, collapse = "_-IN-_"),"_WITHOUT_",(paste0(OutList, collapse= "_-AND-_")))]] <- intersectGR}
      else{holdList[[paste0(paste0(InList, collapse = "_-IN-_"))]] <- intersectGR}
    }
  }
  return(as(holdList, "GRangesList"))  
}
