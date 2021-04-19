ResiliencyGeneticDistance <- function(ReadsTable, RDoutputname, GDoutputname) {

  ReadsTable <- read.csv(ReadsTable, row.names = 1)
  SampleNames <<- colnames(ReadsTable)
 
  TableOfComparisons <- data.frame(vector1 = rep(SampleNames,length(SampleNames)), vector2 = rep(SampleNames, each = length(SampleNames)))

  getrd <- function(location)  {
    vec1name <- location[1]
    vec2name <- location[2]
    print(paste(as.character(vec1name), as.character(vec2name)))
    vec1 <- ReadsTable[,which(SampleNames == vec1name)]
    vec2 <- ReadsTable[,which(SampleNames == vec2name)]
    
    getGD <- function(f1, f2) {
    f1 <- f1/sum(f1)
    f2 <- f2/sum(f2)
    cosTheta=sum(sqrt(f1*f2))
    chorddistance=2*sqrt(2)/pi*sqrt(1-cosTheta)
    chorddistance
    }
    
    times <- min( sum(vec1>0) , sum(vec2>0) )
    
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
    c(rd, gd)
  }
  
  RDtable <- apply(TableOfComparisons, 1, getrd)
  RDtable <- t(RDtable)  
  
  
  RDvector <- as.numeric(RDtable[,1])
  GDvector <- as.numeric(RDtable[,2])
  
  dim(RDvector) <- c(length(SampleNames), length(SampleNames))
  dim(GDvector) <- c(length(SampleNames), length(SampleNames))
  
  RDvector <- as.data.frame(RDvector)
  GDvector <- as.data.frame(GDvector)
  
  rownames(RDvector) <- SampleNames 
  colnames(RDvector) <- SampleNames
  
  rownames(GDvector) <- SampleNames
  colnames(GDvector) <- SampleNames
  
  write.csv(RDvector, RDoutputname)
  write.csv(GDvector, GDoutputname)
  
  RDtable <<- RDtable
  RDvector <<- RDvector
  GDvector <<- GDvector
  ReadsTable <<- ReadsTable
}  
  
  
  
  
  
  
####With Plotting
ResilientGD <- function(vector1, vector2) {
  
 
    vec1name <- vector1
    vec2name <- vector2
    print(paste(as.character(vec1name), as.character(vec2name)))
    vec1 <- ReadsTable[,which(SampleNames == vec1name)]
    vec2 <- ReadsTable[,which(SampleNames == vec2name)]
    
    getGD <- function(f1, f2) {
      f1 <- f1/sum(f1)
      f2 <- f2/sum(f2)
      cosTheta=sum(sqrt(f1*f2))
      chorddistance=2*sqrt(2)/pi*sqrt(1-cosTheta)
      chorddistance
    }
    
    times <- min( sum(vec1>0) , sum(vec2>0) )
    
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
  plot(g, xlab = "Iteration", ylab = "Genetic Distance", ylim = c(0,1))
  print(paste("RD =", rd))
  print(paste("GD =", gd))
  print(paste("MinGD=", mingd))
  print(paste(vec1name, "to", vec2name, "=", log(rd)/log(sum(vec2!=0))))
  print(paste(vec2name, "to", vec1name, "=", log(rd)/log(sum(vec1!=0))))
}
