getFP <- function(ReadsTableName, CFUtable, WhereAreReferences, minweight, outputfilename) {
  library(EnvStats)
  ReadsTable <- read.csv(ReadsTableName, row.names = 1)

  
  
  if(!is.null(CFUtable)) {
    CFUtable <- read.csv(CFUtable)
  }
  
  ReferenceVector <- rowSums(cbind(ReadsTable[,WhereAreReferences]))

 
  
  TestNames <- colnames(ReadsTable)
  TableWithoutNoise <- data.frame(row.names = rownames(ReadsTable))
  
  
  
  ResiliencyIndices <- function(samplename, plots = FALSE, ResiliencyLimit = 20000, FractionSD = 20){
    
    ###Defining Values###
    
    invec <- as.numeric(ReferenceVector)
    outvec <- as.numeric(ReadsTable[,which(TestNames == samplename)])
    
    TotalReads <- sum(outvec)
    outvecfreq <- outvec/sum(outvec)
    outvec[which(outvecfreq < 1E-7)] <- 0
    
    if(!is.null(CFUtable)) {
      cfu <<- CFUtable[CFUtable[,1]==samplename,2]
    } else {cfu <- 1E20}
    
    cfu <- as.numeric(cfu)
    if(length(cfu) == 0) {cfu <- 1E20}

    ###One Resampling Noise Adjustment###
    # ExpectedOnes <- sum(rmultinom(1, sum(outvec), outvec) == 1)
    # ExpectedNonOnes <- sum(rmultinom(1, sum(outvec), outvec) > 1)
    # ExpectedOneRatio <- ExpectedOnes / ExpectedNonOnes
    # 
    # ActualOneRatio <- sum(outvec == 1) / sum(outvec > 1)
    # Correction <- ActualOneRatio / ExpectedOneRatio
    # outvec[outvec <= outvec[outvec > 0][Correction]] <- 0
    # 
    
    ###Actual One ratio Noise Adjustment###
    oneratio <- sqrt( sum(outvec > 1) / sum(outvec == 1))
    if(oneratio < .8) {oneratio <- - (oneratio ^2)}
    NoiseCorrectFactor <- 10^-oneratio
    NoiseDist <- rowSums(ReadsTable[,-WhereAreReferences])
  
    ResampledInvec <- rmultinom(50, round(sum(outvec) * NoiseCorrectFactor), NoiseDist)
    ResampledInvec[ResampledInvec > .005*sum(outvec)] <- 0
    outvecmultiple <- rep(outvec,50)
    dim(outvecmultiple) <- dim(ResampledInvec)
    outvec <- ceiling(rowMeans(outvec - ResampledInvec))
    outvec[outvec < 0] <- 0
    
    ###Defining Times###
   # if ( sum(outvec > 1) >  1.5*sum(outvec == 1) )
   # {times <- min(ceiling(cfu), sum(outvec > 0))
   # } else  {times <- min(ceiling(cfu), sum(outvec > 1))}

    ##smoother between 1 and 1.5

    ratio <- sum(outvec > 1) / sum(outvec == 1)
    factor <- max(0, 1.5-ratio)
    removeones <- min(sum(outvec == 1), round(sum(outvec == 1) * factor))

    if ( sum(outvec > 1) >  1.5*sum(outvec == 1) )
    {times <- min(ceiling(cfu), sum(outvec > 0))
    } else  {times <- min(ceiling(cfu), sum(outvec > 0) - removeones)}

    
    ###Defining greatest frequency change###
    greatestdif <- NULL
  
    outveccopy <- outvec
    outveccopy[outveccopy < 2] <- 0
    outvecsorted <- as.numeric(sort(outveccopy, decreasing = TRUE))
    plotdif <- c(outvecsorted[2:length(outvecsorted)], tail(outvecsorted, 1))
    plotsub <- na.omit(log(outvecsorted) - log(plotdif))
    plotsub <- plotsub[1:length(plotsub)-1]
    
    ###Second Noise Correction###
    if(sum(outvec == 1) > 1.2* sum(outvec > 1)) {
      outvec <- outvec - 1
      greatestdif <- which(plotsub==max(plotsub))
      outvec[outvec < 0] <- 0
      }

    if(max(plotsub) > 2.302585) {greatestdif <- which(plotsub == max(plotsub))}
    
    ###Adjusting times if there are very few nonzero barcodes###
    if (sum(outvec!=0) < 0.1*sum(invec !=0)) {times <- min(ceiling(cfu), sum(outvec > 1))}
      
    ###Starting Resiliency Plot###
    bind <- cbind(invec, outvec)
    bindsorted <- bind[order(bind[,2]),]
    x <- 0
    outvecsorted <- sort(outvec, decreasing = TRUE)
    
    minusone <- function(dims) {
      input <- bindsorted[,1][1:dims]
      out <- bindsorted[,2][1:dims]
      inputprop <- na.omit(input/sum(input))
      outprop <- na.omit(out / sum(out))
      num <- ( outprop - inputprop ) ^ 2
      den <- inputprop*(1-inputprop)
      sigma <- num / den
      sigma <- sigma[sigma!=Inf]
      F <- mean(na.omit(sigma))
      max( 1 / (F - 1/sum(invec) - 1/sum(outvec)), 1/F)
    }
    
    timestorep <- rev(seq(from = length(invec) - times + 1, to = length(invec), by = 1))
    
    if(times > ResiliencyLimit) {x <- rep(1, times)} else 
    {x <- sapply(timestorep, minusone)}
    
    
    scanformin <- function(p) {
      start <- p
      findmin <- function() {
        where <- which(x == start)
        #newlocation <- abs(rnorm(1, mean = where, sd = length(na.omit(x))/10))
        newlocation <- abs(rnorm(1, mean = where, sd = sum(outvec>1)/FractionSD))
        if(newlocation > length(na.omit(x))) {newlocation <- where}
        if(newlocation < 1) {newlocation <- where}
        newstart <- x[round(newlocation)]
        if(newstart < start) {start <- newstart}
        start<<-start
      }
      startvector <- replicate(10000, findmin())
      decision <- which(x==start)
    }
    
    q <- seq(1, length(na.omit(x)), length.out = log(length(na.omit(x))))
    p<-x[q]
    
    if(mean(na.omit(x)) == 1) {guesses <- length(na.omit(x))} else
    {guesses <- sapply(p, scanformin)}
    
    if(length(guesses) == 0) {guesses <- 1}

    ###Defining breaks###
    guessesuniquesorted <- sort(unique(c(guesses)))
    guessesuniquesorted <- c(guessesuniquesorted, length(na.omit(x)))
    
    secondlastguess <- guessesuniquesorted[length(guessesuniquesorted) -1]
    if (rev(sort(outvec))[secondlastguess] - rev(sort(outvec))[max(guessesuniquesorted)] < 4)  {guessesuniquesorted <- guessesuniquesorted[1:length(guessesuniquesorted)-1]}
    
    
     xdif <- c(log(x[1]), log(x))
     xdif <- xdif[1:length(xdif)-1]
     xsub <- (log(x) - xdif)
     greatestdifx <- ifelse(max(xsub) > 2.302585, which(xsub == max(xsub))-1, max(guessesuniquesorted)) 
    
 
    if(is.null(greatestdif)) {greatestdif <- length(na.omit(x))}
    guessesuniquesorted <- sort(unique(c(guessesuniquesorted, greatestdif, greatestdifx)))
    guessesuniquesorted <- guessesuniquesorted[guessesuniquesorted!=0]
    guessesuniquesortedstaggered <- c(-10000000, guessesuniquesorted)
    difference <- c(guessesuniquesorted, 10000000) - guessesuniquesortedstaggered
   
    
     if (min(difference) == 1 & cfu > 2) {
      breaktoremove <- max(which(difference == 1))-1
      guessesuniquesorted <- guessesuniquesorted[-breaktoremove]
    }
    guessesuniquesorted[guessesuniquesorted > max(guessesuniquesorted)-.01*length(na.omit(x))] <- max(guessesuniquesorted)
    guessesuniquesorted <- unique(guessesuniquesorted)
    guessesuniquesorted <<- guessesuniquesorted
    
    ###Creating Indices Table###
    accountsfor <- function(t) {
      outvecsorted <- sort(outvec)
      topnumbers <- tail(outvecsorted, t)
      sum(topnumbers) / sum(outvec)
    }
    
    
    fractionaccounted <- sapply(guessesuniquesorted, accountsfor)
    staggered <- c(0, fractionaccounted)[1:length(fractionaccounted)]
    subtracted <- fractionaccounted - staggered 
    subtracted <- subtracted/sum(subtracted)

    indices <- data.frame(guessesuniquesorted, subtracted)
    colnames(indices) <- c("Number of barcodes", "Fraction of reads")
    
    
    weights <- log(indices[,2])
    values <- (indices[,1])
    weightsforsubtraction <- c(weights[2:length(values)], 0)
    weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
    if(length(weightsdif)==0) {weightsdif <- values}
    noisestart <- indices[,1][which(weightsdif == max(weightsdif))]
    noisestartcopy <- noisestart
    
    #if there are no populations that are below the minimum weight threshold, set noise to be the end of the last detected population
    if (min(indices[,2]) > minweight) {noisestart <- max(indices[,1])}
    #if the sum of the weights of the population after the start of noise is greater than the miniumum weight threshold, set the start to be where the the rest of the reads after are under the minimum weight
    locationofminweightcutoff <- min(which(cumsum(indices[,2])>(1-minweight)))
    if (sum(indices[,2][which(indices[,1] > noisestart)]) > minweight) {noisestart <- indices[,1][locationofminweightcutoff]}
    ##Removing least abundant barcode if it is too nonabundant compared to the second least abundant barcode
    if (noisestart == 2) {
    LeastAbundantBarcode <- sort(outvec, decreasing=TRUE)[noisestart]
    SecondLeastAbundantBarcode <- sort(outvec, decreasing=TRUE)[noisestart-1]
    if(LeastAbundantBarcode/SecondLeastAbundantBarcode < .01) {noisestart <- noisestart - 1}
    }
    
    numberofzeros <- length(invec) - noisestart
    outvecwithoutnoise <- outvec
    outvecwithoutnoise[order(outvecwithoutnoise)][1:numberofzeros] <- 0
    
    ##Removing least abundant barcode if it is too nonabundant compared to the second least abundant barcode
    as.numeric(head(sort(outvecwithoutnoise, decreasing = TRUE)))
    
    FirstResample <- as.numeric(rmultinom(1, sum(outvecwithoutnoise), ReferenceVector/sum(ReferenceVector)))
    
    GetNewBotTable <- function(n) {
      #vec <- as.numeric(rmvhyper(1, FirstResample, n))
      vec <- rmultinom(5, n, FirstResample/sum(FirstResample))
      vec[vec!=0] <- 1
      mean(colSums(vec))
    }
    
    steps1 <- round(seq(from = 1, to = sum(invec != 0), length.out = 100))
    steps2 <- round(seq(from = sum(invec != 0), to = sum(invec != 0)*20, length.out = 100))
    steps <- c(steps1, steps2)
    y <- as.numeric(sapply(steps, GetNewBotTable))

    
    interpol <- approx(x = steps, y = y, n = length(invec))
    xvals <- as.numeric(unlist(interpol[1]))
    yvals <- as.numeric(unlist(interpol[2]))
    
    dfxy2 <- floor(sortedXyData(x = xvals, y = yvals))
    Ns <- NLSstClosestX(dfxy2, noisestart)
    
    inputprop <- na.omit(ReferenceVector/sum(ReferenceVector))
    outprop <- na.omit(outvec / sum(outvec))
    num <- ( outprop - inputprop ) ^ 2
    den <- inputprop*(1-inputprop)
    sigma <- num / den
    sigma <- sigma[sigma!=Inf]
    F <- mean(na.omit(sigma))
    Nb <- 1 / (F - 1/sum(invec) - 1/sum(outvec)) 
    AverageFreq <- 1/geoMean(outvecwithoutnoise[outvecwithoutnoise>0]/sum(outvecwithoutnoise[outvecwithoutnoise>0]))
    
    
    ReadsAtWeightCutoff <- (1-minweight)*sum(outvec)
    MinCutoff <- outvecsorted[min(which(cumsum(outvecsorted[outvecsorted > 0]) > ReadsAtWeightCutoff))]
    NBarcodesAtMinweight <- sum(outvec > MinCutoff) + 1  
    Ns_MinCutoff <- NLSstClosestX(dfxy2, NBarcodesAtMinweight)
    
    if(plots){
      
      converttocutoffs<- function(p) {
        outvecsorted <- sort(outvec, decreasing = TRUE)
        cutoffvalue <- outvecsorted[p]
        if(cutoffvalue == 0) {cutoffvalue <- 1}    
        cutoffvalue/sum(outvec) }
      cutoffpositions <- unlist(lapply(guessesuniquesorted, converttocutoffs))
      cutoffpositions <<- cutoffpositions
      
      par(mfrow = c(2,1))
      plot(1:times,x, log = "y", ylim = c(.01, 2E6), xlim = c(1, length(as.numeric(na.omit(x)))), main = "Resiliency", ylab = "Nb", xlab = "Iteration")
      abline(v=guessesuniquesorted)
      plot(outvec/sum(outvec), log = "y", ylim = c(1E-8, 1), main = samplename, ylab = "Frequency", xlab = "Barcode")
      abline(h = cutoffpositions)
      print(indices)
      print(paste("Nb=", Nb, sep = ""))
      print(paste("Ns=", Ns, sep = ""))
      print(paste("AverageFreq=", AverageFreq, sep = ""))
      print(paste("Ns at minweight=", Ns_MinCutoff, sep= ""))
      
      Nb <<- Nb
      Ns <<- Ns
      AverageFreq <<- AverageFreq
      indices <<- indices
      outvec <<- outvec
      invec <<- invec
      dfxy2 <<- dfxy2
      FirstResample <<- FirstResample
      noisestart <<- noisestart
      x<<-x
      oneratio <<- oneratio
      greatestdif <<- greatestdif
      cfu <<- cfu
      times <<- times
      originaloutvec <<- ReadsTable[,which(TestNames == samplename)]
      MinCutoff <<- MinCutoff
      guesses <<- guesses
    }
    
    
    if(plots == FALSE) {
      print(samplename)
      print(indices)
      print(Ns)
      
      numberofzeros <- length(invec) - noisestart
      outvec[order(outvec)][1:numberofzeros] <- 0
      
      TableWithoutNoise <- data.frame(TableWithoutNoise, outvec)
      colnames(TableWithoutNoise)[length(colnames(TableWithoutNoise))] <- samplename
      TableWithoutNoise <<- TableWithoutNoise
      outvec <<- outvec
      invec <<- ReferenceVector
      if(cfu == 1E20) {cfu <- 0}
      
      den <- density(dfxy2$y)
      denapprox <- approx(den, n=length(invec))
      deny <- denapprox$y
      scalefac <- 1/max(deny)
      uncertainty <- (deny*scalefac)[noisestart]
      
      c(TotalReads, noisestart, Ns_MinCutoff, Nb, Ns, AverageFreq, cfu, log10(Ns), log10(cfu), cfu/Ns)
    }
  }
  

  SampleNames <- colnames(ReadsTable)[-WhereAreReferences]
  TableOfEstimates <- t(sapply(SampleNames, ResiliencyIndices))
  
  colnames(TableOfEstimates) <- c("TotalReads", "Number of barcodes", "Ns_MinCutoff", "Nb", "Ns","AverageFrequency", "CFU", "Log10Ns", "Log10CFU", "CFU/Ns")
  
  write.csv(TableOfEstimates, outputfilename)
  write.csv(TableWithoutNoise, "FrequenciesWithoutNoise.csv")
  TableOfEstimates <<- TableOfEstimates
  FrequenciesWithoutNoise <<- TableWithoutNoise
  ReadsTable <<- ReadsTable
  CFUtable <<- CFUtable
  ReferenceVector <<- ReferenceVector
  TestNames <<- TestNames
  minweight <<- minweight
  ReadsTableName <<- ReadsTableName
  WhereAreReferences <<- WhereAreReferences
}
