
combinereplicates <- function(listoffiles, pvalcutoff, fccutoff, outputfilename) {

combinetolist <- function(filename) {
Finaltable <- read.csv(filename, row.names = 1)
}
combinedlist <- lapply(listoffiles, combinetolist)

combined <- unlist(combinedlist)
numrows <- dim(as.data.frame(combinedlist[1]))[1]
dim(combined) <- c(numrows, length(combined)/numrows)
combined <- as.data.frame(combined)

infcol <- seq(2, 100, 4)[1:length(listoffiles)]
pvalcol <- seq(3, 100, 4)[1:length(listoffiles)]
fccol <- seq(4, 100, 4)[1:length(listoffiles)]

infonly <- combined[,infcol]
pvalonly <- combined[,pvalcol]
fconly <- combined[,fccol]

library(EnvStats)

getlog10mean <- function (x) {
x <- as.numeric(x)
mean(log10(x))
}


meetpvalcutoff <- function(x) {
  x <- as.numeric(x)
  sum(x < pvalcutoff)
}

getlog10sd <- function (x) {
  x <- as.numeric(x)
  sd(log10(x))
}

meanpval <- apply(pvalonly, 1, getlog10mean)
sdpval <- apply(pvalonly, 1, getlog10sd)
passpval <- apply(pvalonly, 1, meetpvalcutoff)
pvalrank <- (length(meanpval) - rank(meanpval))/(length(meanpval)-1)

getsd <- function(x) {
  x <- as.numeric(x)
  sd(x)
}
getmean <- function(x) {
  x<- as.numeric(x)
  mean(x)
}
meetfccutoff <- function(x) {
  x <- as.numeric(x)
  sum(x < fccutoff) + sum(x > abs(fccutoff))
}

meanfc <- apply(fconly, 1, getmean)
sdfc <- apply(fconly, 1, getsd)
passfc <- apply(fconly, 1, meetfccutoff)
fcrank <- (length(meanfc) - rank(-abs(meanfc)))/(length(meanfc)-1)

meaninfsites <- apply(infonly, 1, getmean)

score <- passfc+passpval
score <- score/max(na.omit(score))
score[is.na(score)] <- min(na.omit(score))

rankvector <- rowMeans(cbind(fcrank, pvalrank))
rankvectorfinal <- (length(rankvector)+1) - rank(rowMeans(cbind(rankvector, score)))

master <<- data.frame(combined[,1], meaninfsites, meanpval, sdpval, passpval, meanfc, sdfc, passfc, score, rankvectorfinal)
colnames(master) <<- c("locus tag", "mean inf sites", "mean log10pval", "sd log10pval", paste("pval<", pvalcutoff, sep = ""), "mean log2fc", "sd log2fc", paste("log2fc<", fccutoff, sep=""), "score", "rank")

write.csv(master, outputfilename)

}
