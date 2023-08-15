Updated 10/21/22

This repository contains the code to calculate STAMPR measurements from barcode frequency counts. If you have a question or are having trouble running the code, don't hesitate to email me at hullahalli[at]g.harvard.edu. 

Papers that have deposited code into this repository include:

1) [The description of STAMPR](https://journals.asm.org/doi/10.1128/mSystems.00887-21)
2) [Population dynamics of E. coli systemic infection](https://elifesciences.org/articles/70910)
3) [Barcoded V. cholerae challenge strains](https://journals.asm.org/doi/10.1128/mbio.00539-22)
4) [C. rodentium colonization bottlenecks](https://www.biorxiv.org/content/10.1101/2022.10.11.511778v2)
5) [TLR4 and dose response during E. coli systemic infection](https://www.biorxiv.org/content/10.1101/2023.06.09.543079v1)
-----

The **DoseScalingTLR4** folder contains scripts used [here](https://www.biorxiv.org/content/10.1101/2023.06.09.543079v1)

-----
The **Campbell2023** folder contains the scripts used [here](https://www.biorxiv.org/content/10.1101/2022.10.11.511778v2).

-----
The **ZChol** folder contains a script that is an updated version of the original *getNrNb* function. This was used [here](https://www.biorxiv.org/content/10.1101/2021.12.17.473008v1) to calculate FP values. The original *getNrNb* function runs out of memory when using a large number of barcodes (~60,000 in this study), so this script was created to be able to run on libraries of this size and also improves on noise correction. We anticipate combining the features of both scripts in the future. 

-----

The **STAMPR_Scripts** folder contains all scripts needed to run the analysis from the mSystems STAMPR paper. Note that newer versions have since been deposited in their respective folders

The basic table you need (or will be generated) should be called *ReadsTable*. This is a data frame with the first column specifying the barcode identifier and other columns representing counts for each sample. See *1_to_54_OrderedFrequencies.csv* for a formatting example.


**CompareMutants.R** accepts tables of reads counts across barcodes and outputs the change in barcode frequencies from inocula. 

**ImportSamples.R** takes separate read count files (csv files of one column containing the barcode ID and the second column containing the read count) and merges them into one. Optional plotting function is below the importing function

**IdentifyTopBarcodes.R** is used to compare the top barcodes between two samples. Have one table of reads (output of ImportSamples, plus any reference samples you need to add) called ReadsTable and specify the name of the columns in quotes 

**MajorityDistance.R** computes GD and RD from *ReadsTable* files

**FinalResiliencyWithCFU.R** is the Resiliency Algorithm. Set plots = FALSE if you want to just get Nr/Nb for a bunch of samples, or to TRUE if you want to see the underlying Resiliency plots for a single sample. If you want the latter, set plots to TRUE and load the function called "ResiliencyIndices". All you then need to do is execute the code above the start of the function, and you can run the ResiliencyIndices on any sample in ReadsTable by specifying the column name in quotes. 

-----

The **ExPEC** folder contains the ReadsTable and CFUtable files to reproduce much of our analysis. Some new organ codes are l (entire liver), s (entire spleen) and u (lungs)

*1_to_54_OrderedFrequencies* is the barcode counts for the 32 mice across 12 organs, described in Figure S2-S4 and S6-S11. *CompleteTimecourse_CFU* is the CFU values for these animals. The animal codes are described in the text

*Clod1and2_Frequencies* is barcode counts from the clodonate experiments, and the CFU values are provided as well. 1-4 and 13-16 were given control liposomes, and 6-8 and 17-20 were given clodronate liposomes. 

*Dose3_Frequencies* is barcodes counts from the dosing experiment. 1-10 were given the low dose, and 11-19 were given the high dose. 

*Mut1and2_Frequencies* is barcode counts from Figure 7. 1-8 were from 1 dpi and 9-16 were from 5 dpi.

-----

The **TIS_TAtally** folder contains the read counts across TA sites for all animals in Figure S12 and the replated input. These can be directly input into *CompareTA_NewSim.R* in the **RTISAN_Scripts** folder for comparative analyses. 

