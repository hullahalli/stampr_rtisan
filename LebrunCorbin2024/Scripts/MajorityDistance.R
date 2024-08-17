ResiliencyGeneticDistance <- function(ReadsTableName, expname, ComparisonMetaDataName, limit) {

  ReadsTable <- read.csv(ReadsTableName, row.names = 1)
  MetaData <- read.csv(ComparisonMetaDataName)
  if( length(which(MetaData$Sample == "")) != 0) {print("STOP!!!")}
  
  SampleNames <- MetaData$Sample[order(MetaData$Group, MetaData$Order )]
  
  
  TableOfComparisons <- data.frame(vector1 = rep(SampleNames,length(SampleNames)), vector2 = rep(SampleNames, each = length(SampleNames)))

  FilterToGroups <- function(vec) {
    group1 <- MetaData$Group[which(MetaData$Sample == as.character(vec[1]))]
    group2 <- MetaData$Group[which(MetaData$Sample == as.character(vec[2]))]
    c(group1, group2)
  }
  
  GroupTable <- t(apply(TableOfComparisons,1, FilterToGroups))
  
  FilteredComparisons <- TableOfComparisons[ GroupTable[,1] == GroupTable[,2], ]
  
  getrd <- function(location)  {
    vec1name <- location[1]
    vec2name <- location[2]
    print(paste(as.character(vec1name), as.character(vec2name)))
    vec1 <- ReadsTable[,which(colnames(ReadsTable) == vec1name)]
    vec2 <- ReadsTable[,which(colnames(ReadsTable) == vec2name)]
    
    getGD <- function(f1, f2) {
    f1 <- f1/sum(f1)
    f2 <- f2/sum(f2)
    cosTheta=sum(sqrt(f1*f2))
    if(1-cosTheta < 0) {cosTheta <- 1}
    chorddistance=2*sqrt(2)/pi*sqrt(1-cosTheta)
    chorddistance
    }
    
    times <- min( sum(vec1>0) , sum(vec2>0), limit)
    
    squareroot <- sqrt((vec1/sum(vec1))*(vec2/sum(vec2)))   
    bind <- as.data.frame(cbind(vec1, vec2, squareroot))
    bindsorted <- bind[order(bind[,3]),]
    bindsortedcopy <<- bindsorted
    bindsortedcopy2 <- bindsortedcopy
    
   
    g <- 0
    minusone <- function() {
      vec1 <- bindsorted[,1]
      vec2 <- bindsorted[,2]
      gd <- getGD(vec1, vec2) 
      g <<- c(g, gd)
      bindsorted <<- bindsorted[1:dim(bindsorted)[1]-1,]
    }
    
    replicate(times, minusone())
    g <- g[2:length(g)]
    
    bindsorted <<- bindsorted
    rd <<- sum(na.omit(g) < 0.8) 
    gd <<- g[1]
    c(rd, gd, vec1name, vec2name)
  }
  
  RDtable <- apply(FilteredComparisons, 1, getrd)
  RDtable <- data.frame(t(RDtable))  
  
  
  RDvector <- data.frame(matrix(ncol = length(SampleNames), nrow = length(SampleNames)), row.names = SampleNames)
  colnames(RDvector) <- SampleNames
  RDvector[is.na(RDvector)] = 0
  
  AddToRDMatrix <- function(vec) {
    RDvector[ which (rownames(RDvector) == as.character(vec[3])) , which(rownames(RDvector) == as.character(vec[4])) ] <- vec[1] 
    RDvector <<- RDvector
  }
  
  invisible(apply(RDtable, 1, AddToRDMatrix))
  
  
  GDvector <- data.frame(matrix(ncol = length(SampleNames), nrow = length(SampleNames)), row.names = SampleNames)
  colnames(GDvector) <- SampleNames
  GDvector[is.na(GDvector)] = 0
  
  AddToGDMatrix <- function(vec) {
    GDvector[ which (rownames(GDvector) == as.character(vec[3])) , which(rownames(GDvector) == as.character(vec[4])) ] <- vec[2] 
    GDvector <<- GDvector
  }
  
  invisible(apply(RDtable, 1, AddToGDMatrix))
  
  
  RDvector <- as.data.frame(RDvector)
  GDvector <- as.data.frame(GDvector)
  
  rownames(RDvector) <- SampleNames 
  colnames(RDvector) <- SampleNames
  
  rownames(GDvector) <- SampleNames
  colnames(GDvector) <- SampleNames
  
  ##Correcting RD values for limited iterations
  
  CorrectRD <- function(name) {
    rdvec <- as.numeric(RDvector[,which(colnames(RDvector) == name)])
    position <- which(as.numeric(rdvec) == limit)
    if(sum(position) == 0) {rdvec <- as.numeric(rdvec)}
    else
    {
      NBarcodes <- sum(ReadsTable[,which(colnames(ReadsTable) == name)] != 0)
      rdvec[rdvec == limit] <- NBarcodes
    }
    as.numeric(rdvec)
  }
  
  CorrectedRDvector <- sapply(colnames(RDvector), CorrectRD)
  rownames(CorrectedRDvector) <- SampleNames 
  
  
  ConvertToFRD <- function(vec) {
    vec <- as.numeric(na.omit(vec))
    log(vec+1) / log(max(vec)+1)
  }
  
  CorrectedFRDvector <- apply(CorrectedRDvector, 2, ConvertToFRD)
  
  rownames(CorrectedFRDvector) <- SampleNames
  
  

  write.csv(CorrectedRDvector, paste("CorrectedRD_", expname, ".csv", sep = ""))
  write.csv(GDvector, paste("GD_", expname, ".csv", sep = ""))
  write.csv(CorrectedFRDvector, paste("CorrectedFRD_", expname, ".csv", sep = ""))
  
  RDtable <<- RDtable
  RDvector <<- RDvector
  CorrectedRDvector <<- CorrectedRDvector
  GDvector <<- GDvector
  CorrectedFRDvector <<- CorrectedFRDvector
  ReadsTable <<- ReadsTable
  SampleNames <<- SampleNames
}  
  
  
  
  
  
  
####With Plotting
ResilientGD <- function(vector1, vector2, limit) {
  
 
    vec1name <- vector1
    vec2name <- vector2
    print(paste(as.character(vec1name), as.character(vec2name)))
    vec1 <- ReadsTable[,which(SampleNames == vec1name)]
    vec2 <- ReadsTable[,which(SampleNames == vec2name)]
    
    getGD <- function(f1, f2) {
      f1 <- f1/sum(f1)
      f2 <- f2/sum(f2)
      cosTheta=sum(sqrt(f1*f2))
      if(1-cosTheta < 0) {cosTheta <- 1}
      chorddistance=2*sqrt(2)/pi*sqrt(1-cosTheta)
      chorddistance
    }
    
    times <- min( sum(vec1>0) , sum(vec2>0), limit)
    
    squareroot <- sqrt((vec1/sum(vec1))*(vec2/sum(vec2)))   
    bind <- as.data.frame(cbind(vec1, vec2, squareroot))
    bindsorted <- bind[order(bind[,3]),]
    bindsortedcopy <<- bindsorted
    bindsortedcopy2 <- bindsortedcopy
    
    
    g <- 0
    minusone <- function() {
      vec1 <- bindsorted[,1]
      vec2 <- bindsorted[,2]
      gd <- getGD(vec1, vec2) 
      g <<- c(g, gd)
      bindsorted <<- bindsorted[1:dim(bindsorted)[1]-1,]
    }
    
    replicate(times, minusone())
    g <- g[2:length(g)]
  
    bindsorted <<- bindsorted
    rd <<- sum(na.omit(g) < 0.8)
    gd <<- g[1]
    mingd <- min(na.omit(g))
  
  par(mfrow = c(1,1))
  plot(g, xlab = "Iteration", ylab = "Genetic Distance", ylim = c(0,1), main = paste(vec1name, vec2name ))
  print(paste("RD =", rd))
  print(paste("GD =", gd))
  print(paste("MinGD=", mingd))
  print(paste(vec1name, "to", vec2name, "=", log(rd + 1)/log(sum(vec2!=0) + 1)))
  print(paste(vec2name, "to", vec1name, "=", log(rd + 1)/log(sum(vec1!=0) + 1)))
}
