# this script creates pairwise list comparisons from a vector 
# the options for mode are "pair", "trio", "quadruplet"
# how to prevent pairwise comparison of stuff already compared hmmmmmmm
buildXwiseListFromVec <- function(sampleVec, elementsN = 2)
{
  returnList <- list()
  
    index=1
    for (firstItem in sampleVec)
    {
      if(elementsN ==1)
      {
        return(sampleVec)
      }
      else
      {
        sampleVec2 <- sampleVec[!grepl(firstItem, sampleVec)]
        for (secondItem in sampleVec2)
        {
          if(elementsN == 2)
          {
            if(!(any(duplicated(c(firstItem, secondItem)))))
            {
              returnList[[index]] <- c(firstItem, secondItem)
              index = index +1    
            }
          }
          else
          {
            sampleVec3 <- sampleVec2[!grepl(secondItem, sampleVec2)]
            for (thirdItem in sampleVec3)
            {
              if(elementsN == 3)
              {
                if(any(duplicated(c(firstItem, secondItem, thirdItem)))){}else
                {
                  returnList[[index]] <- c(firstItem,secondItem, thirdItem)
                  index = index +1    
                }
              }
              else
              {
                sampleVec4 <- sampleVec3[!grepl(thirdItem, sampleVec3)]
                for (fourthItem in sampleVec4)
                {
                  if(elementsN == 4)
                  {
                    if(any(duplicated(c(firstItem, secondItem, thirdItem, fourthItem)))){}else
                    {
                      returnList[[index]] <- c(firstItem,secondItem, thirdItem, fourthItem)
                      index = index +1    
                    }
                  }
                  else
                  {
                    sampleVec5 <- sampleVec4[!grepl(fourthItem, sampleVec4)]
                    for (fifthItem in sampleVec5)
                    {
                      if(elementsN == 5)
                      {
                        if(any(duplicated(c(firstItem, secondItem, thirdItem, fourthItem, fifthItem)))){}else
                        {
                          returnList[[index]] <- c(firstItem,secondItem, thirdItem, fourthItem, fifthItem)
                          index = index +1    
                        }
                      }
                      else
                      {
                        sampleVec6 <- sampleVec5[!grepl(fifthItem, sampleVec5)]
                        for (sixthItem in sampleVec6)
                        {
                          if(elementsN == 6)
                          {
                            if(any(duplicated(c(firstItem, secondItem, thirdItem, fourthItem, fifthItem, sixthItem)))){}else
                            {
                              returnList[[index]] <- c(firstItem,secondItem, thirdItem, fourthItem, fifthItem, sixthItem)
                              index = index +1    
                            }
                          }
                          else
                          {
                            print (elementsN)
                          }
                          sampleVec6 <- sampleVec6[!grepl(sixthItem, sampleVec6)]
                        }
                      }
                      sampleVec5 <- sampleVec5[!grepl(fifthItem, sampleVec5)]
                    }
                  }
                  sampleVec4 <- sampleVec4[!grepl(fourthItem, sampleVec4)]
                }
              }
              sampleVec3 <- sampleVec3[!grepl(thirdItem, sampleVec3)]
            }
          }
          sampleVec2 <- sampleVec2[!grepl(secondItem, sampleVec2)]
        }
      }
      sampleVec <- sampleVec[!grepl(firstItem, sampleVec)]
    }
  return(returnList)
}
# for debugging size
#(paste0("Making ", factorial(length(sampleVec))/(factorial(elementsN)*factorial(length(sampleVec)-elementsN)),
      #  " combinations from ", length(sampleVec), " elements, taking ", elementsN, " at a time."))
# sampleList = c("IPSC_H3K27me3", "D7_H3K27me3", "IPSC_H3K27ac", "D7_H3K27ac", "IPSC_yH2AX", "D7_yH2AX")
# elementsN=6
# buildXwiseListFromVec(sampleList, 2)
