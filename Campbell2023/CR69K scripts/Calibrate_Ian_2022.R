CalibrateNs <- function(Standards, SamplesToCalibrate, outputfilename) {
  standards <- read.csv(Standards)
  samples <- read.csv(SamplesToCalibrate, row.names = 1)
  

  cfu <- log10(standards[,1])
  ns <- log10(standards[,2])
  
  
  fit <- smooth.spline(x = cfu, y = ns, df = 7)
  xvals <- seq(from = min(cfu), to = 7, by = .01)
  
  spline <- predict(fit, x = xvals)
  x_cfu <- as.numeric(unlist(spline[1]))
  y_ns <- as.numeric(unlist(spline[2]))
  
  dfxy <- sortedXyData(x = x_cfu, y = y_ns)
  max_dfxy_x <- NLSstClosestX(dfxy, max(dfxy$y))
  dfxy$y[dfxy$x > max_dfxy_x] <- max(dfxy$y)
  
  
  cal <- function(k) {
    if(k > min(x_cfu)) {10^(NLSstClosestX(dfxy, k))
    } else {10^k}
      
  }
  
  Ns_calibrated <- sapply(log10(samples$Ns), cal)
  caltable <- cbind(samples, Ns_calibrated)
  write.csv(caltable, outputfilename)
  
  
 #plot(predict(fit, x = xvals), cex = .1, ylim = c(0, 6), xlim = c(0, 8), xlab = "CFU from standards (blue) or calibrated Ns (red)", ylab = "Uncalibrated Ns")
 plot((dfxy), cex = .1, ylim = c(0, 6), xlim = c(0, 8), xlab = "CFU from standards (blue) or calibrated Ns (red)", ylab = "Uncalibrated Ns")
 points(cfu,ns, col = "blue")
 points(y = log10(caltable$Ns), x = log10(caltable$Ns_calibrated), col = "red")
 dfxy <<- dfxy
 
}
