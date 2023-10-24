###For NextSeq 1000, strengthRow = .02 and strengthCol = .01 worked well
###For Miseq, strengthRow = .01 and strengthCol = .005 worked well

GetNoHoppingTable <- function(ReadsTableName, SampleTableName, OutputFileName, strengthRow, strengthCol, Scale) {
  
  ReadsTable <- read.csv(ReadsTableName, row.names = 1)
  
  samplenames <- colnames(ReadsTable)
  samplenames <- samplenames[samplenames != "barcode"]
  
  GenerateRowNoiseDistributions <- function(SampleTableName) {
    SampleTable <- read.csv(SampleTableName, row.names = 1)
    
    GetDist <- function(name) {
      row <- which(SampleTable == name, arr.ind = TRUE)[1]
      col <- which(SampleTable == name, arr.ind = TRUE)[2]
      ColsToCombine <- (1:ncol(SampleTable))[-col]
      CombinedIndexes <- cbind(row, ColsToCombine)
      
      GetName <- function(indexes) {
        SampleTable[indexes[1],indexes[2]]
      }
      
      SamplesToCombine <- apply(CombinedIndexes,1,GetName)
      
      SubMatrix <- ReadsTable[,colnames(ReadsTable) %in% SamplesToCombine]
      if(class(SubMatrix) == "integer") {SubMatrix <- cbind(SubMatrix, rep(0, length(SubMatrix)))}
      rowSums(SubMatrix)
      
    }
    
    NoiseDist <- data.frame(sapply(samplenames, GetDist))
    NoiseDist$barcode <- ReadsTable$barcode
    
    
    colnames(NoiseDist) <- na.omit(colnames(ReadsTable))
    
    
    ReadsTableRowNoiseOnly<<-NoiseDist
  }
  
  GenerateRowNoiseDistributions(SampleTableName)
  
  GenerateColNoiseDistributions <- function(SampleTableName) {
    SampleTable <- read.csv(SampleTableName, row.names = 1)
    
    GetDist <- function(name) {
      row <- which(SampleTable == name, arr.ind = TRUE)[1]
      col <- which(SampleTable == name, arr.ind = TRUE)[2]
      RowsToCombine <- (1:nrow(SampleTable))[-row]
      CombinedIndexes <- cbind(RowsToCombine, col)
      
      
      GetName <- function(indexes) {
        SampleTable[indexes[1],indexes[2]]
      }
      
      SamplesToCombine <- apply(CombinedIndexes,1,GetName)
      
      SubMatrix <- ReadsTable[,colnames(ReadsTable) %in% SamplesToCombine]
      rowSums(SubMatrix)
      
    }
    
    NoiseDist <- data.frame(sapply(samplenames, GetDist))
    NoiseDist$barcode <- ReadsTable$barcode
    
    
    colnames(NoiseDist) <- na.omit(colnames(ReadsTable))
    
    
    ReadsTableColNoiseOnly<<-NoiseDist
  }
  
  
  GenerateColNoiseDistributions(SampleTableName)
  
 if ( sum(ReadsTableRowNoiseOnly == 0) == dim(ReadsTableRowNoiseOnly)[1] * dim(ReadsTableRowNoiseOnly)[2]) {ReadsTableRowNoiseOnly[ReadsTableRowNoiseOnly == 0] <- .0001}

  
  SubtractNoise <- function(strengthRow, strengthCol, OutputFileName){
    
    ##First Two scale with SAMPLE
    
    if(Scale == "sample") {
    AdjustForStrengthRow <- function(name) {
      vecoriginal <- ReadsTable[,which(colnames(ReadsTable) == name)]
      vec <- ReadsTableRowNoiseOnly[,which(colnames(ReadsTableRowNoiseOnly) == name)]
      strength <- round(sum(vecoriginal)) * strengthRow
      factor <- strength / sum(vec)
      vec * factor
    }

    AdjustForStrengthCol <- function(name) {
      vecoriginal <- ReadsTable[,which(colnames(ReadsTable) == name)]
      vec <- ReadsTableColNoiseOnly[,which(colnames(ReadsTableColNoiseOnly) == name)]
      strength <- round(sum(vecoriginal)) * strengthCol
      factor <- strength / sum(vec)
      vec * factor
    }
    }
    
    #These two scale with NOISE
    
    if(Scale == "noise"){
    AdjustForStrengthRow <- function(name) {
      vec <- ReadsTableRowNoiseOnly[,which(colnames(ReadsTableRowNoiseOnly) == name)]
      strength <- round(sum(vec)) * strengthRow
      factor <- strength / sum(vec)
      vec * factor
    }
    
    AdjustForStrengthCol <- function(name) {
      vec <- ReadsTableColNoiseOnly[,which(colnames(ReadsTableColNoiseOnly) == name)]
      strength <- round(sum(vec)) * strengthCol
      factor <- strength / sum(vec)
      vec * factor
    }
    }
    
    NoiseAdjustedRow <- data.frame(sapply(samplenames, AdjustForStrengthRow))
    NoiseAdjustedCol <- data.frame(sapply(samplenames, AdjustForStrengthCol))
    
    colnames(NoiseAdjustedRow) <- samplenames
    colnames(NoiseAdjustedCol) <- samplenames
    
    NoiseAdjusted <- NoiseAdjustedRow + NoiseAdjustedCol
    
    RemoveNoise <- function(name) {
      vecoriginal <- ReadsTable[,which(colnames(ReadsTable) == name)]
      vecnoise <- NoiseAdjusted[,which(colnames(NoiseAdjusted) == name)]
      ans <- round(vecoriginal - vecnoise)
      ans[ans < 0] <- 0
      ans
    }
    
    ReadsTableNoHopping <- data.frame(sapply(samplenames, RemoveNoise))
    ReadsTableNoHopping <<- data.frame("barcode" = rownames(ReadsTable), ReadsTableNoHopping)
  }
  
  SubtractNoise(strengthRow = strengthRow, strengthCol=strengthCol, OutputFileName = OutputFileName)
  
  
  
  GetNoiseReport <- function(name) {
    sum(ReadsTableNoHopping[,which(colnames(ReadsTableNoHopping) ==name )]) /  sum(ReadsTable[,which(colnames(ReadsTable) ==name )]) 
    
  }
  
  NoiseReport <- data.frame(sapply(samplenames, GetNoiseReport))
  colnames(NoiseReport)
  
  ReadsTable <<- ReadsTable
  ReadsTableNoHopping <<- ReadsTableNoHopping
  samplenames <<- samplenames
  strengthRow <<- strengthRow
  strengthCol <<- strengthCol
  Scale <<- Scale
  write.csv(ReadsTableNoHopping, OutputFileName, row.names = FALSE)
  write.csv(NoiseReport, "NoiseReport.csv")
}

plotSTAMPNoHopping <- function(name) {
  tiffname <- paste(name, "_nohopping.tiff", sep="")
  tiff(tiffname, width = 1000, height = 1000, points = 20)
  column <- which(name == colnames(ReadsTableNoHopping))
  x <- 1:dim(ReadsTableNoHopping)[1]
  y <- log10(ReadsTableNoHopping[,column]/sum(ReadsTableNoHopping[,column]))
  plot(x, y, ylim = c(-8, 0), xlab = "Barcode", ylab = "Log10 Frequency", main = paste(name, "NoHopping", " RowStrength = ", strengthRow, " ColStrength = ", strengthCol, Scale, sep = " "))
  dev.off()}

sapply(samplenames, plotSTAMPNoHopping)






compareplotsNoise <- function(nameofref, samplesize) {
  
  reflocationRow <- which(colnames(ReadsTableRowNoiseOnly) == nameofref)
  reflocationCol <- which(colnames(ReadsTableColNoiseOnly) == nameofref)
  
  lookuplocation <- which(colnames(ReadsTable) == nameofref)
  
  referenceRow <- ReadsTableRowNoiseOnly[,reflocationRow]
  referenceCol <- ReadsTableColNoiseOnly[,reflocationCol]
  lookup <- ReadsTable[,lookuplocation]
  
  referencepropRow <- referenceRow/sum(referenceRow)
  referencepropCol<- referenceCol/sum(referenceCol)
  lookupprop  <- lookup/sum(lookup)
  
  combined <- as.data.frame(cbind(referencepropRow, referencepropCol, lookupprop))
  
  combinedsortRow <- combined[order(combined[,1]),]
  toprefvaluesRow <- tail(combinedsortRow, samplesize)
  
  combinedsortCol <- combined[order(combined[,2]),]
  toprefvaluesCol <- tail(combinedsortCol, samplesize)
  
  
  par(mfrow = c(3, 1))
  
  plot(log10(referenceRow/sum(referenceRow)), yaxt = "none", ylim = c(-6,0), ylab = "Log10 Frequency", xlab = "Barcode", cex = .5, main = paste(nameofref, "Row Noise"))
  
  points(as.numeric(rownames(toprefvaluesRow)), log10(toprefvaluesRow[,1]), cex = 1, col = "red", pch = 3)
  
  points(as.numeric(rownames(toprefvaluesCol)), log10(toprefvaluesCol[,1]), cex = 1, col = "blue", pch = 4)
  
  axis(2, seq(-6,0,1),las=2)
  
  
  plot(log10(referenceCol/sum(referenceCol)), yaxt = "none", ylim = c(-6,0), ylab = "Log10 Frequency", xlab = "Barcode", cex = .5, main = paste(nameofref, "Col Noise"))
  
  points(as.numeric(rownames(toprefvaluesRow)), log10(toprefvaluesRow[,2]), cex = 1, col = "red", pch = 3)
  
  points(as.numeric(rownames(toprefvaluesCol)), log10(toprefvaluesCol[,2]), cex = 1, col = "blue", pch = 4)
  
  axis(2, seq(-6,0,1),las=2)
  
  
  
  
  plot(log10(lookup/sum(lookup)), yaxt = "none", ylim = c(-6,0), ylab = "Log10 Frequency", xlab = "Barcode", cex = .5, main = nameofref)
  
  points(as.numeric(rownames(toprefvaluesRow)), log10(toprefvaluesRow[,3]), cex = 1, col = "red", pch = 3)
  
  points(as.numeric(rownames(toprefvaluesCol)), log10(toprefvaluesCol[,3]), cex = 1, col = "blue", pch = 4)
  
  axis(2, seq(-6,0,1),las=2)
  
  
  
}
