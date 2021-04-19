# stampr_rtisan

Updated 4/19/2021

Maintenance of this repository will be continuous as our manuscripts are in revision. Files will be updated and more will be added as needed. If you are here from our bioRxiv manuscripts and have a question, feel free to email me at hullahalli[at]g[dot]harvard[dot]edu. 


The STAMPR_Scripts folder contains all scripts needed to run the analysis. 

CompareMutants.R accepts tables of reads counts across barcodes and outputs the change in barcode frequencies from inocula. 

ImportSamples.R takes separate read count files (csv files of one column containing the barcode ID and the second column containing the read count) and merges them into one. Optional plotting function is below the importing function

IdentifyTopBarcodes.R is used to compare the top barcodes between two samples. Have one table of reads (output of ImportSamples, plus any reference samples you need to add) called ReadsTable and specify the name of the columns in quotes 

MajorityDistance.R computes GD and RD from ReadsTable files

FinalResiliencyWithCFU.R is the Resiliency Algorithm. Set plots = FALSE if you want to just get Nr/Nb for a bunch of samples, or to TRUE if you want to see the underlying Resiliency plots for a single sample. 





The ExPEC folder contains the ReadsTable and CFUtable files to reproduce much of our analysis. We set InocCFU = 2E7, minweight = .03, CorrectForNoise = .005. Feel free to vary these parameters to see how it changes the results. References are columns 1-5
