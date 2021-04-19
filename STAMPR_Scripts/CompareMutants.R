CompareMutants <- function(ReadsTable, WhereAreInocula, outputfilename) {
  
  ReadsTable <- read.csv(ReadsTable, row.names = 1)
  RefTable <- ReadsTable[,WhereAreInocula]
  SampleTable <- ReadsTable[,-WhereAreInocula]
  
  RefVector <- rowSums(RefTable)
  RefVectorProp <- RefVector / sum(RefVector)
  
  DivideWithinColumns <- function(q) {
    q[q==0] <- 1
    CI <- (q/sum(q))  / RefVectorProp
  }
  
  NormalizedCI <- apply(SampleTable, 2, DivideWithinColumns)
 
  write.csv(NormalizedCI, outputfilename)
  
  
} 
### trying to transform this but its not worth the effort
TranformTable <- function(columns) {
    Cols <- grepl(columns, colnames(NormalizedCI))
    SubsetOrgan <- as.matrix(NormalizedCI[,Cols])
    
    AllMutants <- matrix()
    
    GetRows <- function(rows) {
      r <- grepl(rows, rownames(SubsetOrgan))
      SubsetMutant <- as.matrix(SubsetOrgan[r,])
      
      namesofrows <- rownames(SubsetMutant)
      namesofcols <- colnames(SubsetMutant)  
      product <- length(namesofrows)*length(namesofcols)
    
      rownamerep <- rep(namesofrows, length.out = product)
      colnamerep <- rep(namesofcols, 1, each = length(namesofrows))
    
      bind <- cbind(rownamerep, colnamerep)
      
      combinenames <- function(p) {
        paste(p[1], p[2], sep = "-")  
      }
      
      newnames <- apply(bind, 1, combinenames)  
      
      dim(SubsetMutant) <- c(1, length(SubsetMutant))
      
      NewSubsetMutant <- data.frame(newnames, as.numeric(SubsetMutant))
      
      colnames(NewSubsetMutant) <- c("Pair", rows)
      
      
      AllMutants <- data.frame(AllMutants, NewSubsetMutant)
      
      AllMutants <<- AllMutants
      
      
    }
    
   lapply(mutants, GetRows)
   
   
    
    
    
    
  }
  
  
  
  
  write.csv(NormalizedCI, outputfilename)
  
}
