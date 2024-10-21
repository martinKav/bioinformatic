# list of small functions
filterByScoreGR <- function(GR, quant){ GR[score(GR) >= quantile(score(GR), probs=quant)]}
filterByTotalSignal <- function(GR, quant){ GR[GR$TotalSignal >= quantile(GR$TotalSignal, probs=quant)]}
mcolsRemove <- function(columnName){if(grepl("mcols", columnName)){newName <- strsplit(columnName, "mcols.")[[1]][2]; return(newName)}else{return(columnName)}}