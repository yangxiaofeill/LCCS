########################################################
## Local cluster of Cytosine (LCCS) detection
########################################################
library(data.table)
rm(list = ls())

args = commandArgs(trailingOnly=TRUE)
printHelpMessage <- function(){
  cat(paste("#################################################################################### \n", 
            "This program is used to detect Local Clustering of Cytosine Sites (LCCS) of \n", 
            "single nucletide resolution methylation sequencing data \n",
            "Usage: Rscript --vanilla LCCS.R minCNum maxSize posFile \n",
            "Parameters: \n",
            "minCNum, the minimal Cytosine Number in a LCCS, default value is 3 \n",
            "maxSize, the maxmize size of a LCCS, default value is 500bp \n",
            "posFile, the file includes positions of all cytosine, and at lest two column, like \n" ,
            "chr1 1 \n",
            "chr1 5 \n",
            "... \n",
            "#####################################################################################\n",
            "Reference: doi: 10.1093/hmg/ddv172 ",
            "copyright @ Xiaofei Yang, xfyang@xjtu.edu.cn", sep = '')
  )
  stop("", call.=FALSE)
}

printEndMessage <- function(theBedFile){
  newBedFile <- paste(theBedFile, ".LCCS.bed", sep = '')
  newBedFile2 <- paste(theBedFile, ".isolated_cytosine.bed", sep = '')
  cat(paste("#################################################################################### \n", 
            "the prgram finished, and generated LCCS file in a bed format (", newBedFile, "), like: \n",
            "chr1 1 100 LCCS_1 5 \n",
            "the first column is chromosomes, the second column is the start position (0-based), \n ", 
            "the thrid column is the end position (1-based), the forth column is the LCCS ID, \n", 
            "and the fifth column is the number of Cytosines in this LCCS\n",
            "---------------------------------------------------------------------------------- \n",
            "also generated a file included all Cytosines that do not belong to any LCCS (", newBedFile2, "), in the input bed format \n",
            sep = '')
  )
}

theCytosinePosFile <- ""
minCNum <- 3
maxSize <- 500

if(length(args) == 0){
  printHelpMessage()
}else if(length(args) == 1){
  theCytosinePosFile <- args[1]
}else if(length(args) == 3){
  minCNum <- as.numeric(args[1])
  maxSize <- as.numeric(args[2])
  theCytosinePosFile <- args[3]
}else{
  printHelpMessage()
}

LCCS <- function(leftIndex, rightIndex, theData, currLCCS, maxSize, isInLCCS){
  leftPos <- theData$V2[leftIndex]
  rightPos <- theData$V2[rightIndex]
  theSize <- rightPos - leftPos
  
  leftUpIndex <- leftIndex - 1
  leftUpPos   <- theData$V2[leftUpIndex]
  rightDownIndex <- rightIndex + 1
  rightDownPos   <- theData$V2[rightDownIndex]
  
  leftDis <- leftPos - leftUpPos
  rightDis <- rightDownPos - rightPos
  
  if(leftIndex == 1){
    leftDis <- maxSize
    leftUpIndex <- leftIndex
  }
  if(isInLCCS[leftUpIndex] == 1){
    leftDis <- maxSize
  }
  if(rightIndex == nrow(theData)){
    rightDis <- maxSize
    rightDownIndex <- rightIndex
  }
  if(isInLCCS[rightDownIndex] == 1){
    rightDis <- maxSize
  }
  
  if(leftDis <= rightDis & (theSize + leftDis) <= maxSize){
    leftIndex <- leftUpIndex
    currLCCS <- c(currLCCS, leftIndex)
    currLCCS <- LCCS(leftIndex, rightIndex, theData, currLCCS, maxSize, isInLCCS)
  }else if(rightDis < leftDis & (theSize + rightDis) <= maxSize){
    rightIndex <- rightDownIndex
    currLCCS <- c(currLCCS, rightIndex)
    currLCCS <- LCCS(leftIndex, rightIndex, theData, currLCCS, maxSize, isInLCCS)
  }else{
    ## an LCCS generated
    return(currLCCS)
  }
}

getLCCS <- function(chr, theDistance, theData, theCytosinePosFile){
  ALLLCCS <- data.frame(chr = '', start = 0, end = 0, id = "", C_num = 0)
  isInLCCS <- rep(0, nrow(theData))
  theDistanceOrder <- order(theDistance)
  LCCS_ID <- 0
  for(oneOrder in theDistanceOrder){
    currLCCS <- c()
    oneDis <- theDistance[oneOrder]
    if(oneDis <= maxSize){
      leftIndex <- oneOrder
      rightIndex <- oneOrder + 1
      
      isOk <- FALSE
      if(isInLCCS[leftIndex] == 0){
        currLCCS <- c(currLCCS, leftIndex)
        isOk <- TRUE
      }
      if(isInLCCS[rightIndex] == 0){
        currLCCS <- c(currLCCS, rightIndex)
        isOk <- TRUE
      }
      
      if(isOk){
        currLCCS <- LCCS(leftIndex, rightIndex, theData, currLCCS, maxSize, isInLCCS)
        cNum <- length(currLCCS)
        if(length(currLCCS) >= minCNum){
          currLCCS <- sort(currLCCS)
          isInLCCS[currLCCS] <- 1
          LCCS_ID <- LCCS_ID + 1
          ALLLCCS <- rbind(ALLLCCS, data.frame(chr = chr, start = (theData$V2[currLCCS[1]] - 1), 
                                               end = (theData$V2[currLCCS[cNum]]), 
                                               id = paste("LCCS_", LCCS_ID, sep = ''), 
                                               C_num = cNum))
        }
      }
    }
  }
  ALLLCCS <- ALLLCCS[-1, ]
  isolatedC <- theData[which(isInLCCS == 0), ]
  return(list(lccs = ALLLCCS, cyto = isolatedC))
}

print("reading bed file...")
theCPos <- fread(theCytosinePosFile, header = F, stringsAsFactors = F)
theChrs <- unique(theCPos$V1)
ALLLCCS <- data.frame(chr = '', start = 0, end = 0, id = "", C_num = 0)
isolatedC <- theCPos[1, ]
print("detecing the LCCS ...")
for(oneChr in theChrs){
  oneChrData <- theCPos[which(theCPos$V1 == oneChr), ]
  oneChrData <- oneChrData[order(oneChrData$V2), ]
  CNum <- nrow(oneChrData)
  theDistance <- oneChrData$V2[-1] - oneChrData$V2[-CNum]
  oneChrRes <- getLCCS(chr = oneChr, theDistance = theDistance, theData = oneChrData)
  ALLLCCS <- rbind(ALLLCCS, oneChrRes[["lccs"]])
  isolatedC <- rbind(isolatedC, oneChrRes[['cyto']])
}
ALLLCCS <- ALLLCCS[-1, ]
isolatedC <- isolatedC[-1, ]
fwrite(ALLLCCS, file = paste(theCytosinePosFile, ".LCCS.bed", sep = ''), sep = "\t", row.names = F, col.names = F, quote = F)
fwrite(isolatedC, file = paste(theCytosinePosFile, ".isolated_cytosine.bed", sep = ''), sep = "\t", row.names = F, col.names = F, quote = F)

printEndMessage(theCytosinePosFile)
