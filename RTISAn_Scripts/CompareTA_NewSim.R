CompareTA_NewSim<- function(TAtally_in, TAtally_out, filename, graphs) {
  
#Imports files
  TAtally_in <- read.csv(TAtally_in)
  TAtally_out <- read.csv(TAtally_out)
  
#Read normalization. Script is run twice to work on both inputs and outputs
print("normalizing input")
library(evobiR)

#This script only normalizes using the bottom 99th percentile, which ignores the effect of genes with unusually high TA hits.
mean99 <- function(x) {
  top <- as.numeric(quantile(x, .99))
  x <- x[x<top]
  mean(x)
}

#Identifies means across each sliding window, set default to 6000. Remember this is 6000 columns of the TAtally table, not the genome. This corresponds to approximately 6000 * 16 = 96000bp of genome
slidingmean <- SlidingWindow(mean99, TAtally_in$Reads, 6000,6000)
slidingmeannozero <- slidingmean[slidingmean != 0]
slidingmean[slidingmean == 0] <- min(slidingmeannozero)

#Identifies the minimum value of each sliding window mean, and corrects each mean to the minimum. This adjustment value is replicated so it matches the length of the vector of TAs, and then simply multiplied. 
slidingmeanadj <- min(slidingmean) / slidingmean
factor <-rep(slidingmeanadj, each = 6000, length.out = length(TAtally_in$Reads))
TAnormin <<- factor*TAtally_in$Reads

print("normalizing output")
slidingmean <- SlidingWindow(mean99, TAtally_out$Reads, 6000,6000)
slidingmeannozero <- slidingmean[slidingmean != 0]
slidingmean[slidingmean == 0] <- min(slidingmeannozero)
slidingmeanadj <- min(slidingmean) / slidingmean
factor <-rep(slidingmeanadj, each = 6000, length.out = length(TAtally_out$Reads))
TAnormout <<- factor*TAtally_out$Reads
if(is.na(sum(TAnormout == TRUE))) {TAnormout <- TAtally_out$Reads}

#Simulating input
print("identifying simulation sampling size")

#First part identifies how many sites to sample so that fraction of TA sites hit in the output equals the fraction of TA sites hit in the simulation
library(extraDistr)

#Assigns variables for splitting the TAnorm_input into 100 segments and the sum of all reads in the input. It's rounded because the hypergeometric distribution doesn't like decimals
z <- round(sum(round(TAnormin)) / 100)
a <- sum(round(TAnormin))

#Generates 100 numbers evenly spaced out starting from 1000 anding with the total number of reads in TAnorm_input
t <- as.list(seq(1000, a, z))

#Creates a function that identifies fraction of TAsites hit after sampling from a multivariate hypergeometric distrubution
getgeomarray <- function (p) {
  sum(rmvhyper(1, round(TAnormin), p) != 0) / length(TAnormin)
}

#Performs the previously made function on the "t" vector made before. The output becomes the Y axis, which tells you the fraction of TA sites hit as a function of how many times you samplied the distrbution
points<- unlist(lapply(t, getgeomarray))

#Identifies the point on the X axis (number of times to sample) that equals the fraction of TAsites hit in the outpu
dfxy <- sortedXyData(unlist(t), points)
colnames(dfxy) <- c("x", "y")
n <- sum(TAnormout != 0) / length(TAnormout)
readstosample <- NLSstClosestX(dfxy, n)

print("refining estimate")
#This function helps refine the estimate of the number of reads to sample. Basically, it takes the guess made before, and then hops around that value with a normal distrbution, comparing the value it gets to the original value. If the value of the new sampling number is closer to what we want, it replaces the other number, and if not, it tries again. This goes on and on and we take the mean and sd sampling number of 200 tries. These are then used to get a normal distribution from which we sample a vector used to do the actual simulation
refineestimate <- function() {
  
  startingvalue <- as.vector(rmvhyper(1, round(TAnormin), round(readstosample)))
  startingfrac <- sum(startingvalue != 0) / length(startingvalue)
  #These next 4 lines try to guess a good standard deviation to use for the sampling, and I've right now set it equal to the sd of the local region of the dfxy table where you are trying to bottleneck to times the actual reads to sample. For example, if your initial guess for readstosample from before was 30,000, this part finds which fraction TAs hit this corresponds to from dfxy, and makes an array from the 5 values above and below it. The sd of the new vector is equal to the SD of the array of 10 values (which are fractionTAsites hit) times 30,000. This means that a sample with a wider bottleneck will have a smaller sd, because the curve of sampling size vs fractionTAsites hit is relatively flat for high sampling sizes (which equal wide bottlenecks). 
  start <- max(which(dfxy$x < readstosample))
  up <- start+5
  down<- abs(start-5)
  range <- as.numeric(na.omit(dfxy$y[down:up]))
  newsampling <- rnorm(1, mean = readstosample, sd = sd(range*readstosample))
  if (newsampling > sum(round(TAnormin))) {newsampling <- sum(round(TAnormin))}

  #To see if the new function is a good guess, we take the average fractionTAsites hit from 5 runs of the hypergeometric 
  newvaluematrix <- t(rmvhyper(5, round(TAnormin), round(newsampling)))
  newfracvector <- colSums(newvaluematrix != 0) / length(TAnormin)
  newfrac <<- mean(newfracvector)
  
  #the readstosample variable is changed only if the next guess is closer to the desired value than the previous one
  if(abs(newfrac-n) < abs(startingfrac-n)) {readstosample <- newsampling} 
  if(readstosample > sum(round(TAnormin))) {readstosample <- sum(round(TAnormin))}
  readstosample <<- readstosample
}

#This is repeated 200 times, and then then we can get a normal distribution of sampling reads
refinedvector <<- replicate(200, refineestimate())
readstosamplemean <- mean(refinedvector)
readstosamplesd <- sd(refinedvector)
arraytosample <- rnorm(100, mean = readstosamplemean, sd = readstosamplesd)
arraytosample[arraytosample > sum(round(TAnormin))] <- sum(round(TAnormin))

#Performs the simulation 100 times, and finally multiplicatively scales so that the number of reads between each column of the simulation and the output are equal.
print("simulating input")

simfinal <- function(x) {
  as.numeric(rmvhyper(1, round(TAnormin), round(x)))
}
sim <- lapply(arraytosample, simfinal)
sim <- unlist(sim)
dim(sim) <- c(length(TAnormin), 100)
simsums <- colSums(sim)
simcorrect <- sum(TAnormout)/simsums
siminput <<- as.data.frame(sweep(sim, MARGIN = 2, simcorrect, "*"))


#Begin Analysis

#Combines the TA position, reads from output, and simulated input into one matrix for simplicity
  master <<- cbind(TAtally_in$Location, TAnormout, siminput)
  colnames(master) <<- c("Location", "TAnormout", c(1:100))

#Removes all intergenic regions from the master matrix, including TA sites for which the mean of all values in all columns is 0 - these correspond to TAs that have no insertions in the input or the output  
  masterNoIg <- master[master$Location != "IG",]
  masterNoIgNoZero <- masterNoIg[which(rowMeans(masterNoIg[,c(2:102)])!=0),]

#Generates vector and character strings corresponding the output reads and the locus tags. There will of course be many redundancies in the "Locations" string, since there are many TA sites per gene
  TAnormoutNoIg <- as.vector(masterNoIgNoZero$TAnormout)
  Locations <- as.character(masterNoIgNoZero$Location)

#Extracts unique locus tags and separates the simulation from the new master matrix that lacks IGs and TA sites with 0 reads
  LocustagsUnique<-unique(masterNoIgNoZero$Location)
  masterSimNoIG <- unname(as.matrix(masterNoIgNoZero[,3:102]))
  
#Begin Mann Whitney
  print("performing mann-whitney tests")
  library(matrixTests)
  library(EnvStats)

#This function performs the Mann-Whitney. It looks for where x, which in the downstream lapply function will be a locus tag, is located, first in the "Locations" vector of all nonzero, nonIG TA sites. Using these position, it extracts the corresponding values in the TAnorm_output and corresponding rows in nonzero nonIG simulation. These two are inputs for a column by column fast MWU test. NAs are then removed. It is critical that the length of "TAnormoutNoIG" equals the number of rows of "masterSimNoIG". The if statement here is a techniality of the col_wilcoxon_twosample function - for any gene that is one row long(only one nonzero TAsite), we add an additional row so that the function actually runs column by column, since it only does that on samples with multiple rows. These will not pass significance anyway.
  wilcoxtest <- function(x) {
    spots <- as.vector(which(Locations == x))
    a<-masterSimNoIG[spots,]
    z<-TAnormoutNoIg[spots]
    if (class(a) == "numeric") {
      a<-rbind(a,a) 
      z<-rbind(z,z)
      }
    h <- col_wilcoxon_twosample(a,z)$pvalue
    h[!is.na(h)]
  }

#This lapply function goes through the list of unique nonzero locus tags and performs the function above. The results get stored as a list, which then gets unlisted after comupting the geometric mean of the P-values
  MWUtable<-lapply(LocustagsUnique, wilcoxtest)
  MWUmeans <- unlist(lapply(MWUtable, geoMean))
 
#A curation step to correct the artifical P values - basically sets them to the max p value.   
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  modePval <- Mode(MWUmeans)
  MWUmeans[MWUmeans == modePval] <- max(MWUmeans)

  print(paste("P-val 0.05 quantile =", as.numeric(quantile(MWUmeans, 0.05))))

  print("calculating fold changes") 
#This function sums all the number of reads at each TA site in each column of the simulation for a given gene, specified as "x" and is input from lapply below. The output is a vector that contains the sum of reads at a given gene for all runs of the simulation
  simSums <- function(x) {
    spots <- as.vector(which(Locations == x))
    a<-masterSimNoIG[spots,]
    if(class(a) == "matrix") { colSums(a) } else { a }
  }

#Lapply to do this across all unique nonzero locus tags, and then a second lapply to average the sums across all simulations.
  Atable <- lapply(LocustagsUnique, simSums)
  Ameans <- unlist(lapply(Atable, mean)) 
  
#This is vey similar to the previous step, except it's performed on the TAnorm_output vector. Because it's just a vector, it's much simpler
  outputSums <- function(x) {
    spots <- as.vector(which(Locations == x))
    sum(TAnormoutNoIg[spots])
  }
  
  Ztable <- unlist(lapply(LocustagsUnique, outputSums))
  
#Replaces zeros in the output file, which would otherwise give values of 0 for a fold change, with the minimum detectable value observed in the table of sums for the output.
  minval <- min(Ztable[Ztable!=0])
  Ztable[Ztable==0] <- minval

#Divides the sums of the TAnorm_output file with the mean of the sums of the simulation to calculate a fold change.
  FCtable <- Ztable / Ameans
  
  print(paste("Log2FC 0.05 quantile =", log2(as.numeric(quantile(FCtable, 0.05)))))
  
  print("calculating informative sites")
#The informative site calculation requires a different set of variables, since we do want to keep TA sites that have no reads.  Intergenic regions are still discarded. We also set every nonzero value to one, so that the sum is equal to the number of informative sites
  sitearray <- (masterNoIg[,2:102])
  sitearray[sitearray != 0] <- 1
  Locustagssitearray <- masterNoIg$Location
  Locustagsinf <- unique(masterNoIg$Location)

#Similar to before, we set x to be a locus tag and look for x's position across our locus tag array. This is then crossed with the matrix of 1s and 0s at each TAsite. The colums are summed to indicate how many informative sites are present in the output and 100 simulations. Then we take the average of this.
  infsites <- function(x) {
    spots <- as.vector(which(Locustagssitearray == x))
    mean(as.vector(colSums(sitearray[spots,])))
  }
  
#lapply is performed just like before.
  infarray <- unlist(lapply(Locustagsinf, infsites))
  
#Resulting values of FC, MWU, and informative sites are combined. Since the informative sites vector will be longer (because we kept genes with 0 insertions), we need to use the merge function and impute NAs for where there is no corresponding P-value or fold change.
  Amerge<- data.frame(Locustagsinf, infarray)
  colnames(Amerge) <- c("Locus_tag", "Inf_sites")
  Bmerge <- data.frame(LocustagsUnique, MWUmeans, log(FCtable,2))
  colnames(Bmerge) <- c("Locus_tag", "MWU_P-value", "Log2FC")
  Finaltable <<- merge(Amerge, Bmerge, all.x = TRUE)
  write.csv(Finaltable, filename)
  
#Few other things are output to help with troubleshooting and visualizing your data, which is an optional argument
  TAtally_in <<- TAtally_in
  TAtally_out <<- TAtally_out
  notzero <- function(x) {
    sum(x != 0) / length(TAnormin) }
  simulateddistribution <<- apply(siminput, 2, notzero)
  
  if(graphs) {
  print("generating plots")
  par(mfrow = c(3,2))
  
  uniquelocations <- unique(TAtally_in$Location)
  uniquelocations <- uniquelocations[uniquelocations != "IG"]
  
  fracthitinput <- function(x) {
    subset <- TAtally_in[which(TAtally_in$Location == x),]
    sum(subset$Reads !=0 ) / length(subset$Reads)
  }
  
  fracinput <<- unlist(lapply(uniquelocations, fracthitinput))
  
  fracthitoutput <- function(x) {
    subset <- TAtally_out[which(TAtally_out$Location == x),]
    sum(subset$Reads !=0 ) / length(subset$Reads)
  }
  
  fracoutput <<- unlist(lapply(uniquelocations, fracthitoutput))
  
  hist(fracoutput, col=rgb(0,0,1,0.5), breaks = 50, main = "input(red) vs output(blue) hit fractions", ylab = "Frequency", xlab = "Fraction TA sites hit")
  hist(fracinput, col=rgb(1,0,0,0.5),  breaks = 50, add = T)
  
  
  hist(simulateddistribution, main = "Simulated complexities", xlab = "Fraction TA sites hit", breaks = 10)
  
  plot(x = TAtally_out$TA/1E6, y = TAtally_out$Reads, log = "y", ylab  = "number of reads", xlab = "TA position (Mb)", main = "Original output Reads")
  plot(x = TAtally_out$TA/1E6, y = TAtally_in$Reads, log = "y", ylab  = "number of reads", xlab = "TA position (Mb)", main = "Original input Reads")
  plot(x = TAtally_out$TA/1E6, y = TAnormout, log = "y", ylab  = "number of reads", xlab = "TA position (Mb)", main = "Normalized Output Reads")
  plot(x = TAtally_out$TA/1E6, y = TAnormin, log = "y", ylab  = "number of reads", xlab = "TA position (Mb)", main = "Normalized input Reads")
  }
  
}
