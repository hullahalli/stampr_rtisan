getNrNb <- function(ReadsTable, CFUtable, InocCFU, WhereAreReferences, minweight, CorrectForNoise, outputfilename, CalibrationFile = NULL) {
  ##Ensure that plots is set to false - this toggle then lets you run the resiliency script alone to look at things in more fine resolution
  
  #These four lines set up the necessary metadata, and split up the table of reads into inputs and outputs 
  ReadsTable <- read.csv(ReadsTable, row.names = 1)
  
  if(!is.null(CFUtable)) {
  CFUtable <- read.csv(CFUtable)
  }
  
  ReferenceVector <- rowMeans(cbind(ReadsTable[,WhereAreReferences]))
  TestNames <- colnames(ReadsTable)[-WhereAreReferences]
  TableWithoutNoise <- data.frame(row.names = rownames(ReadsTable))
  
 #Makes a table for estimating bottleneck from fraction of barcodes that are identified
  
  library(extraDistr)
  RoundedRefVector <-unname(round(round(ReferenceVector) * InocCFU / sum(round(ReferenceVector))))
  steps <- seq(from = 1, to = length(RoundedRefVector)*10, by = 10)
  
  GetBotTable <- function(n) {
    vec <- as.numeric(rmvhyper(1, RoundedRefVector, n))
    sum(vec!=0)
  }
  
  y <- unlist(lapply(steps, GetBotTable))
  
  dfxy <- sortedXyData(steps, y)
  colnames(dfxy) <- c("x", "y")
  
  #Internal Resiliency Function - this is the bulk of the script. The function gets applied on the vector of sample names (TestNames above), so this is run for every noninput sample in your dataset
  ResiliencyIndices <- function(samplename, plots = TRUE){
    
    #Specifies input vector (which is the average of your inputs), output (the row which corresponds to the particular sample name), and CFU (a single value corresponding to the name of the sample)
    invec <- ReferenceVector
    outvec <- ReadsTable[,colnames(ReadsTable) == samplename]
    
    if(!is.null(CFUtable)) {
    cfu <<- CFUtable[CFUtable[,1]==samplename,2]
    } else {cfu <- 1E8}
    
    #Noise adjustment
    if (CorrectForNoise != 0) {
      ResampledInvec <- rmvhyper(1, round(invec), round(CorrectForNoise*sum(outvec)))
      outvec <- as.numeric(outvec - ResampledInvec)
      outvec[outvec < 0] <- 0
    }
    
    
    #Times specifies how many iterations the resiliency function is run. Too many or too few isn't ideal, but we can constrain it by CFU. So if you have 2 CFU, it will only run 2 times, because it's pointless to try to look farther than that. If you have a lot of CFU, this variable is set to the number of barcodes that have > 1 read. 1 is arbitrary, but we can be really quite confident that anything 1 or below is noise.
    times <- min(round(cfu+1), sum(outvec > 1))
    
    #Combines input and output vectors into one data frame, and then orders everything by the output vector. The barcodes with the most reads in the output will be at the bottom of the data frame
    bind <- as.data.frame(cbind(invec, outvec))
    bindsorted <- bind[order(bind[,2]),]
    
    #This logical puts in a 0 for all barcodes in the output where we initially define noise. 
    #if(times < length(outvec)) {bindsorted[,2][1:(length(outvec)-times)] <- 0}
    
    #As the script is run, bindsorted itself will be changed, so we make a few copies of it for later use
    bindsortedcopy <<- bindsorted
    bindsortedcopy2 <- bindsortedcopy
    
    #x is the first iteration of the Resiliency graph. The first value is 0 so it has a place to start
    x <- 0
    
    #This function calculates Nb and then appends the value to the x variable. At the end, it shortens the bindsorted script to remove the barcode with the highest read in the output.
    minusone <- function() {
      input <- bindsorted[,1]
      out <- bindsorted[,2]
      inputprop <- na.omit(input/sum(input))
      outprop <- na.omit(out / sum(out))
      num <- ( outprop - inputprop ) ^ 2
      den <- inputprop*(1-inputprop)
      sigma <- num / den
      sigma <- sigma[sigma!=Inf]
      F <- mean(na.omit(sigma))
      nb <- 1 / (F - 1/sum(invec) - 1/sum(outvec)) ##note that this Nb is assuming that the total number of reads as bindsorted gets trimmed is always equal to the total number of reads in the initial sample = otherwise you get negative Nbs and the x variable gets thrown off
      x <<- c(x, nb)
      bindsorted <<- bindsorted[1:dim(bindsorted)[1]-1,]
      
    }
    
    #minusone is run for the number of times specified by times. The initial 0 is removed from x and both x and the new bindsorted are placed into the environment
    replicate(times, minusone())
    x <- x[2:length(x)]
    x<<-x
    bindsorted <<- bindsorted
    
    
    #Outer layer of the function that scans for minima. The input to this, p, is later defined. p is a specific set of numbers that specifiy all the poisitions for which the minimum finder will start. scanformin is applied over each of these positions
    scanformin <- function(p) {
      start <- p
      
      #This is finds minima. It takes p, the initial guess, and looks around with a specified sd (the newlocation variable). It asks if this new value (newstart)  is less than start. If so, start is set to newstart. Repeating this function iteratively changes the value of start if the newstart guess is smaller
      findmin <- function() {
        where <- which(x == start)
        newlocation <- abs(rnorm(1, mean = where, sd = length(na.omit(x))/10))
        if(newlocation > length(na.omit(x))) {newlocation <- where}
        if(newlocation < 1) {newlocation <- where}
        newstart <- x[round(newlocation)]
        if(newstart < start) {start <- newstart}
        start<<-start
      }
      
      #startvector is the resulting output of the findmin function run 1000 times. The start variable also changes to be the last guess. The position of this guess in the x variable is defined as decision. Put another way, decision is the x coordinate of the resiliency graph, and start is the y coordinate.
      startvector <- replicate(1000, findmin())
      decision <- which(x==start)
    }
    
    #Here we define p, which first requires defining q, which is set to be 1/15 of the total number of barcodes. So if this is is 1500 numbers long, q will be 1500/15 = 100 numbers long (1, 15, 30, 45 ... 1500). Correspondingly, p is the value of x at each of these positions. 
    
    q <- round(seq(1, length(na.omit(x)), length.out = length(na.omit(x))/15))
    p<-x[q]
    
    #The guesses variable is defined as all the final decisionsreached from starting across all the values of p.  The unique guesses are taken and sorted.
    guesses <<- unlist(lapply(p, scanformin))
    guessesuniquesorted <- sort(unique(c(guesses)))
    guessesuniquesorted <- c(guessesuniquesorted, length(na.omit(x)))
    
    
    #There is sometimes some weirdness at the end of the run, and it makes decisions really close together. This part takes anything within 5 values of the end of the run and sets it to the end. This new sorted list is uniqued again.
    guessesuniquesorted[guessesuniquesorted > max(guessesuniquesorted)-5] <- max(guessesuniquesorted)
    guessesuniquesorted <- unique(guessesuniquesorted)
    
    
    #Identifies greatest log change
    xdif <- c(log(x[1]), log(x))
    xdif <- xdif[1:length(xdif)-1]
    xsub <- (log(x) - xdif)
    greatestdif = which(xsub == max(xsub)) -1
    guessesuniquesorted <- sort(unique(c(guessesuniquesorted, greatestdif)))
    guessesuniquesorted <- guessesuniquesorted[guessesuniquesorted!=0]
    
    #Sometimes, guesses will be within one unit and this can be problematic when it sees this as a big weight jump (because only one barcode is in that group). This part makes sure that each guess is at least 2 barcodes apart
    
    guessesuniquesortedstaggered <- c(-10000, guessesuniquesorted)
    difference <- c(guessesuniquesorted, 10000) - guessesuniquesortedstaggered
    if (min(difference) == 1) {
      breaktoremove <- max(which(difference == 1)) - 1
      guessesuniquesorted <- guessesuniquesorted[-breaktoremove]
    }
    guessesuniquesorted <<- guessesuniquesorted
    
    
    
    #Once we have the guesses, we need to make the table that describes the weights assigned to each guess. This function takes as an input the guessesuniquesorted list above and identifies the fraction of reads acounded for by all barcodes that have more reads than specific by each value in guessesuniquesorted. For example, if one of the guesses was 15, and position 15 was a barcode with 100 reads, this function will determine the fraction of all reads acounted for by barcodes with at least 100 reads.
    
    accountsfor <- function(t) {
      outvecsorted <- sort(outvec)
      topnumbers <- tail(outvecsorted, t)
      sum(topnumbers) / sum(outvec)
    }
    fractionaccounted <- unlist(lapply(guessesuniquesorted, accountsfor))
    
    #Here we adjust the fractionaccounted variable so that it takes in the fraction of accounted reads in between guesses, not just from that guess upwards    
    staggered <- c(0, fractionaccounted)[1:length(fractionaccounted)]
    subtracted <- fractionaccounted - staggered 
    
    #Here we obtain the max Nb up to each location defined in guessesuniquesorted  
    getintervalnb <- function(y) {
      max(x[1:y])
    }
    nbintervals <- unlist(lapply(guessesuniquesorted, getintervalnb))
    
    #And finally we combine these into the indices table. The headings explain what each value means. We have now subseted our data into descrete segments separated by local minima    
    indices <<- data.frame(nbintervals, subtracted, guessesuniquesorted)
    colnames(indices) <<- c("Nr", "Accounts for", "Number of barcodes")
    
    #These logicals define the start of noise (noisestart) as the location of greatest log change in the "accounts for" section, with a few more housekeeping stuff 
    weights <- log(indices[,2])
    values <- (indices[,1])
    weightsforsubtraction <- c(weights[2:length(values)], 0)
    weightsdif <- (weights - weightsforsubtraction)[1:length(weights)-1]
    if(length(weightsdif)==0) {weightsdif <- values}
    noisestart <<- indices[,3][which(weightsdif == max(weightsdif))]
    noisestartcopy <- noisestart
    
    #if there are no populations that are below the minimum weight threshold, set noise to be the end of the last detected population
    if (is.na(indices[,3][min(which(indices[,2] < minweight))])) {noisestart <<- max(indices[,3])}
    #if the sum of the weights of the population after the start of noise is greater than the miniumum weight threshold, set the start to be where the the rest of the reads after are under the minimum weight
    locationofminweightcutoff <- min(which(cumsum(indices[,2])>(1-minweight)))
    if (sum(indices[,2][which(indices[,3] > noisestart)]) > minweight) {noisestart <<- indices[,3][locationofminweightcutoff]}
    
    
    #Plots x and barcodes if specified   
    if(plots) {
      par(mfrow = c(3,1))
      plot(1:times,x, log = "y", ylim = c(.01, 2E6), xlim = c(1, length(as.numeric(na.omit(x)))), main = "Resiliency", ylab = "Nb", xlab = "Iteration")
      abline(v=guessesuniquesorted)
      plot(outvec/sum(outvec), log = "y", ylim = c(1E-6, 1), main = "Barcodes", ylab = "Frequency", xlab = "Barcode")
    }
    
    
    #This function takes the values in guessesunquesorted and identifies their location in the plot of barcode frequencies vs barcode, which is defined by cutoffpositions
    converttocutoffs<- function(p) {
      cutoffposition <- length(outvec) - p
      cutoffvalue <- bindsortedcopy[cutoffposition+1,2]
      if(cutoffvalue == 0) {cutoffvalue <- 1}    
      cutoffvalue/sum(outvec) }
    
    cutoffpositions <- unlist(lapply(guessesuniquesorted, converttocutoffs))
    cutoffpositions <<- cutoffpositions
    
    #Now we redefine our output vector such that everything after the start of noise is set to 0. Note that this is done on bindsortedcopy2, and not bindsorted or bindsortedcopy
    lengthofnoise <- dim(bindsortedcopy2)[1]-noisestart
    bindsortedcopy2[1:lengthofnoise,][2] <- 0
    outvecwithoutnoise <<- as.numeric(bindsortedcopy2[,2])
    
    #The resiliency function is run again, this time with the new noise-less output vector and where the resulting output is defined as z. z is functionally the same as x, except z is done after the noise correction. 
    z <- 0
    minusonefinal <- function() {
      input <- bindsortedcopy2[,1]
      out <- bindsortedcopy2[,2]
      inputprop <- na.omit(input/sum(input))
      outprop <- na.omit(out / sum(out))
      num <- ( outprop - inputprop ) ^ 2
      den <- inputprop*(1-inputprop)
      sigma <- num / den
      sigma <- sigma[sigma!=Inf]
      F <- mean(na.omit(sigma))
      #
      nb <- 1 / (F- 1/sum(invec) - 1/sum(outvec)) 
      z <<- c(z, nb)
      bindsortedcopy2 <<- bindsortedcopy2[1:dim(bindsortedcopy2)[1]-1,]
    }
    replicate(sum(bindsortedcopy2[2] != 0), minusonefinal())
    z <- z[2:length(z)]
    z<<-z
    bindsortedcopy2 <<- bindsortedcopy2
    
    #Nr is the defined by final value. It is further constrained to ensure that that it is greater than Nb (recall that the original Nb estimate is x[1])   
    finalvalue <<- max(z) 
    if(finalvalue < x[1]) {finalvalue <<- x[1]}
    
    #Adjustment for when Nb < the identified number of barcodes above noise. Here, the number of barcodes is a better estimate. So we change Nr to reflect this.
    minimumnearesttimes <- max(indices[,3][(indices[,3] <= noisestart)])
    CompBottleneck <<- NLSstClosestX(dfxy, minimumnearesttimes)
    if(finalvalue < CompBottleneck) {finalvalue <<- CompBottleneck}
    
    #Plots the new resiliency function if needed
    if(plots){
      abline(h = cutoffpositions)
      print(indices)
      print(finalvalue)
      plot(z, log = "y", ylim = c(.01, 2E6), xlim = c(1, length(as.numeric(na.omit(z)))), main = "Resiliency without noise", ylab = "Nb", xlab = "Iteration")
    }
    
    #If plots is set to FALSE, this indicates to the computer that you are running a lot of data to get our the table, and so this will be gathered.
    if(plots == FALSE) {
      print(samplename)
      if(!is.null(CalibrationFile)) {LookUpTable <<- read.csv(CalibrationFile, row.names = 1)}
      if(x[1] < 1) {x[1] <- 1}
      if(finalvalue<1) {finalvalue <- 1}
      if(!is.null(CalibrationFile)) {
        nbcalibrated <<- as.numeric(LookUpTable[as.character(round(log10(x[1]), digits = 2)),])
        nrcalibrated <<- as.numeric(LookUpTable[as.character(round(log10(finalvalue), digits = 2)),])
      }
      if(is.null(CalibrationFile)) {
        nbcalibrated <- NA
        nrcalibrated <- NA
      }
      
      numberofzeros <- sum(outvecwithoutnoise == 0)
      outvec[order(outvec)][1:numberofzeros] <- 0
      TableWithoutNoise <- data.frame(TableWithoutNoise, outvec)
      colnames(TableWithoutNoise)[length(colnames(TableWithoutNoise))] <- samplename
      TableWithoutNoise <<- TableWithoutNoise
      outvec <<- outvec
      invecNoiseCorrected <<- invec
      invec <<- ReferenceVector
      c(x[1], finalvalue, nbcalibrated, nrcalibrated)
    }
  }
  
  #TableOfEstimates is what the final table is called, and it will have different dimentions depending on if you specify a calibration curve. This variable is written into your directory and important metadata is spit out of the function so you can run the resiliency script to stand alone if needed (remembering to set plots = TRUE)
  TableOfEstimates <- t(sapply(TestNames, ResiliencyIndices))
  if(!is.null(CalibrationFile)) {colnames(TableOfEstimates) <- c("Nb", "Nr", "NbLower_Cal", "NbMedian_Cal", "NbUpper_Cal", "NrLower_Cal", "NrMedian_Cal", "NrUpper_Cal") }
  if(is.null(CalibrationFile)) {colnames(TableOfEstimates) <- c("Nb", "Nr", "Nb_cal", "Nr_cal") }
  
  write.csv(TableOfEstimates, outputfilename)
  write.csv(TableWithoutNoise, "FrequenciesWithoutNoise.csv")
  
  TableOfEstimates <<- TableOfEstimates
  FrequenciesWithoutNoise <<- TableWithoutNoise
  ReadsTable <<- ReadsTable
  CFUtable <<- CFUtable
  ReferenceVector <<- ReferenceVector
  TestNames <<- TestNames
  minweight <<- minweight
  dfxy <<- dfxy
  
}
