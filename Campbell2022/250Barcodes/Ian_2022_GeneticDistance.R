##################################################################################
###                               USER  INPUT                                  ###

### define working directory
setwd("~/Desktop/current_lab_members/Ian/STAMP/V2 20210413/20200915_CR240STAMP_5d_new/add master")

### define input sequence data file ("BigMatAll.csv")
inputdata="master_frequencies.csv"

### read input sequence data file
Data=as.matrix(read.csv(inputdata,row.names=1))

### import calibration data ("CalibrationLookup")
#load("D:/Anleitungen/Methoden Computer/Illumina Sequencing/Bottleneck_analysis_tools/R/R_create_calibration_curve/CalibrationLookup.Rdata")

### define output filename
name="ChordDistance"

### define which columns of BigMatAll contain real read numbers and not metadata like binaries, avarge count number or absolute presence in samples. Usually column 1 to 24 (1:24) contain real read numbers.
v.reads=1:68

### define which columns of BigMatAll contain reference samples (usually the inoculum)? Define at least two consequtive columns e.g. 1:2 would mean column 1 and column 2 of BigMatAll contain reference samples. If there is only one reference sample, repeat the same column (e.g. 1:1)
v.reference=1:4

### USER ### define which columns of BigMatAll contain data samples (that should be the columns of BigMatAll that contain real read numbers minus the reference columns)
v.rabbitsamples=4:68

###                           END OF USER INPUT                                 ###
###################################################################################


### define threshold for inocula (e.g. sequence must be present at least once)
v.threshold=which(rowSums(Data[,v.reference],na.rm=TRUE)>0)

### choose real reads
ThresholdMatrix=Data[v.threshold,v.reads]

### NA= 0 (not present)
ThresholdMatrix[is.na(ThresholdMatrix)==T]=0

### Matrix mit nur Referenzen machen
ReferenceMatrix=ThresholdMatrix[,v.reference]

###FrequenzMatrix fuer Referenzen definieren
fReferenceMatrix=matrix(nrow=nrow(ReferenceMatrix),ncol=ncol(ReferenceMatrix),dimnames=list(rownames(ReferenceMatrix),colnames(ReferenceMatrix)))

###FrequenzMatrix mit Frequenzen fuellen
for(x in 1:ncol(fReferenceMatrix)){
	fReferenceMatrix[,x]=ReferenceMatrix[,x]/sum(ReferenceMatrix[,x])
}

### Average frequencies
fTotalInoculum=rowMeans(fReferenceMatrix,na.rm=T)

### Make Matrix for Frequencies
fRabbitMatrix=matrix(nrow=length(fTotalInoculum),ncol=length(v.rabbitsamples))
colnames(fRabbitMatrix)=c(colnames(Data[,v.rabbitsamples]))

RabbitMatrix=ThresholdMatrix[,v.rabbitsamples]
for(y in 1:ncol(fRabbitMatrix)){
	fRabbitMatrix[,y]=RabbitMatrix[,y]/sum(RabbitMatrix[,y])
}
fTotalMatrix=cbind(fReferenceMatrix,fRabbitMatrix)

### DistanceMatrix
RabbitDistance=matrix(ncol=ncol(fTotalMatrix),nrow=ncol(fTotalMatrix),dimnames=list(colnames(fTotalMatrix),colnames(fTotalMatrix)))
for(i in 1:ncol(fTotalMatrix)){
	for(ii in 1:ncol(fTotalMatrix)){
		f1=fTotalMatrix[,i]
		f2=fTotalMatrix[,ii]
		cosTheta=sum(sqrt(f1*f2))
		theta=acos(cosTheta)/(2*pi)*360
		chorddistance=2*sqrt(2)/pi*sqrt(1-cosTheta)
		RabbitDistance[i,ii]=chorddistance
	}
}

library(ape)
tree=nj(as.dist(RabbitDistance))
treephylo=as.phylo(tree)
plot(tree)

###DistannzMatrix angucken
RabbitDistance
write.csv(RabbitDistance,file="DistanceMatrix.csv")
###image(x=colnames(fRabbitMatrix),y=colnames(fRabbitMatrix),log10(RabbitDistance))

### clears workspace
rm(list=ls())

