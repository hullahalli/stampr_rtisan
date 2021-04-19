compareplots <- function(nameofref, nameoflookup, samplesize) {
  
  
  reflocation <- which(colnames(ReadsTable) == nameofref)
  lookuplocation <- which(colnames(ReadsTable) == nameoflookup)
  
  reference <- ReadsTable[,reflocation]
  lookup <- ReadsTable[,lookuplocation]
  
  referenceprop <- reference/sum(reference)
  lookupprop  <- lookup/sum(lookup)
  combined <- as.data.frame(cbind(referenceprop, lookupprop))
  combinedsort <- combined[order(combined[,1]),]
  toprefvalues <- tail(combinedsort, samplesize)
  
  
  par(mfrow = c(2, 1))
  
  plot(log10(reference/sum(reference)), yaxt = "none", ylim = c(-6,0), ylab = "Log10 Frequency", xlab = "Barcode", cex = .5, main = nameofref)

  points(as.numeric(rownames(toprefvalues)), log10(toprefvalues[,1]), cex = 1, col = "red", pch = 10)
  axis(2, seq(-6,0,1),las=2)
  
  
  plot(log10(lookup/sum(lookup)), yaxt = "none", ylim = c(-6,0), ylab = "Log10 Frequency", xlab = "Barcode", cex = .5, main = nameoflookup)
  
  points(as.numeric(rownames(toprefvalues)), log10(toprefvalues[,2]), cex = 1, col = "red", pch = 10)
  axis(2, seq(-6,0,1),las=2)
  

}
