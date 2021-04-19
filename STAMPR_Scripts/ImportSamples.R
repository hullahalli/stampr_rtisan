importfiles <- function(name) {
  
namesoffiles <- list.files()
importedlist <- lapply(namesoffiles, read.csv)
unlisted <- unlist(importedlist)
numberofsamples <- as.numeric(length(namesoffiles))
numberofrows <- numberofsamples * 2
numberofbarcodes <- as.numeric(length(unlisted) / numberofsamples / 2)
dim(unlisted) <- c(numberofbarcodes, numberofrows)
columnswithreads <- seq(from = 2, to = dim(unlisted)[2], by = 2)
tableofreads <- as.data.frame(unlisted[,columnswithreads])
tableofreads <- apply(tableofreads, 2, as.numeric)
barcodelist <- unlisted[,1]

takename <- function(x) {
  y <- unlist(strsplit(x, " "))
  y[1]
  
}

samplenames <<- unlist(lapply(namesoffiles, takename))
namedtable <- data.frame(barcodelist, tableofreads)
colnames(namedtable) <- c("barcode", samplenames)
BarcodeFrequencies <<- namedtable
write.csv(BarcodeFrequencies, name, row.names = FALSE)
}


plotSTAMP <- function(name) 
 {tiffname <- paste(name, ".tiff", sep="")
  tiff(tiffname, width = 1000, height = 1000, points = 20)
  column <- which(name == colnames(BarcodeFrequencies))
  x <- 1:dim(BarcodeFrequencies)[1]
  y <- log10(BarcodeFrequencies[,column]/sum(BarcodeFrequencies[,column]))
  plot(x, y, ylim = c(-6, 0), xlab = "Barcode", ylab = "Log10 Frequency", main =name)
  dev.off()}

sapply(samplenames, plotSTAMP)





###merging all 5 pools together

firstmaster <- read.csv("FirstPoolFrequencies.csv", row.names = 1)
secondmaster <- read.csv("SecondPoolFrequencies.csv", row.names = 1)
thirdmaster <- read.csv("ThirdPoolFrequencies.csv", row.names = 1)
fourthmaster <- read.csv("FourthPoolFrequencies.csv", row.names = 1)
fifthmaster <- read.csv("FifthPoolFrequencies.csv", row.names = 1)




addcolumns <- function (df1, df2) {
  df1names <- colnames(df1)
  df2names <- colnames(df2)
 
  start <- df1$barcode
  
  addnewcolumns <- function(t) {
    if(length(which(df1names == t)) == 0) {start <- data.frame(start, as.numeric(df2[,df2names ==t]))}
    
    if(length(which(df1names == t)) == 1) {start <- data.frame(start, rowSums(cbind(as.numeric(df2[,df2names ==t]), as.numeric(df1[,df1names ==t]))))}
    colnames(start)[dim(start)[2]] <- t
    start <<- start
  }
  
  sapply(df2names, addnewcolumns)
  
  reinsertcolumns <- function(t) {
    if(length(which(df2names == t)) == 0) {
      start <- data.frame(start, as.numeric(df1[,df1names ==t]))
      colnames(start)[dim(start)[2]] <- t
      }
    
    start <<- start
  }
  
  sapply(df1names, reinsertcolumns)
  
  
  combineddf <<- start
}

addcolumns(firstmaster, secondmaster)
addcolumns(combineddf, thirdmaster)
addcolumns(combineddf, fourthmaster)
addcolumns(combineddf, fifthmaster)

combineddf <- data.frame(firstmaster$barcode, combineddf)  

write.csv(combineddf, "FinalMaster.csv")
##now manually curate and reimported, removed in.p2

