Updated 4/30/2021

Maintenance of this repository will be continuous as our manuscripts are in revision and beyond. Files will be updated and more will be added as needed. If you have a question or are having trouble running the code, feel free to email me at hullahalli[at]g.harvard.edu. 

The scripts in this repository are from two manuscripts. The [first](https://www.biorxiv.org/content/10.1101/2021.04.28.441820v1) is a computational one describing STAMPR, a methodology to improve quantification of bottlenecks and dissemination patterns during infection. The other is a biological one describing E. coli systemic infection population dynamics, which will be uploaded in a few weeks. If you are reading this now, we are sort of in a strange middle ground where one of the manuscripts is available online but the other isn't. 

The animal experiments from the biology manuscript gave us a good set of data to develop STAMPR. To avoid future data duplication problems, the computational manuscript doesn't talk about any systemic ExPEC biology, but simply uses examples from this data to better demonstrate how the algorithm works. *1_to_54_OrderedFrequencies.csv* in the **ExPEC** folder contains this data to facilitate review. Each column of this file represents a different biological sample, but for now, these just serve as a bunch of test samples with which the algorithm can be run. Our goal is to split up these manuscripts in a way where you can read one if you want to learn about STAMPR and the other if you want to learn about ExPEC systemic infection.

-----

The **STAMPR_Scripts** folder contains all scripts needed to run the analysis. 
The basic table you need (or will be generated) should be called *ReadsTable*. This is a data frame with the first column specifying the barcode identifier and other columns representing counts for each sample. See *1_to_54_OrderedFrequencies.csv* for a formatting example.


**CompareMutants.R** accepts tables of reads counts across barcodes and outputs the change in barcode frequencies from inocula. 

**ImportSamples.R** takes separate read count files (csv files of one column containing the barcode ID and the second column containing the read count) and merges them into one. Optional plotting function is below the importing function

**IdentifyTopBarcodes.R** is used to compare the top barcodes between two samples. Have one table of reads (output of ImportSamples, plus any reference samples you need to add) called ReadsTable and specify the name of the columns in quotes 

**MajorityDistance.R** computes GD and RD from *ReadsTable* files

**FinalResiliencyWithCFU.R** is the Resiliency Algorithm. Set plots = FALSE if you want to just get Nr/Nb for a bunch of samples, or to TRUE if you want to see the underlying Resiliency plots for a single sample. If you want the latter, set plots to TRUE and load the function called "ResiliencyIndices". All you then need to do is execute the code above the start of the function, and you can run the ResiliencyIndices on any sample in ReadsTable by specifying the column name in quotes. 

-----

The **ExPEC** folder contains the ReadsTable and CFUtable files to reproduce much of our analysis. We set InocCFU = 2E7, minweight = .03, CorrectForNoise = .005. Feel free to vary these parameters to see how it changes the results. References are columns 1-5
