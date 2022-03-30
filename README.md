# Estimation of serial intervals using viral sequences

This repository contains code for estimating the serial interval distribution in case clusters using genomic data, taking incomplete sampling into account.

This is supporting code for <br />
**'Genomic epidemiology offers high resolution estimates of serial intervals for COVID-19'** <br />
Jessica E. Stockdale, Kurnia Susvitasari, Paul Tupper, Benjamin Sobkowiak, Nicola Mulberry, Anders Gonçalves da Silva, Anne E. Watt, Norelle Sherry, Corinna Minko, Benjamin P. Howden, Courtney R. Lane, and Caroline Colijn.

## File/folder contents:
* 'SI-from-clusters' contains all user-facing code and the R project file. A readme with full details of all scripts is included in this folder.
* 'SI with coprimary' contains utility functions for the underlying statistical model including coprimary transmission.
* 'SI with noncoprimary' contains equivalent utilities for a model with no coprimary transmission.
* 'data_simulated' contains a simulated dataset which is ready to be used with the above code
 
 All code was prepared and tested under R version 4.1.0 (2021-05-18). Package dependencies are listed and will be automatically installed in *SI-from-clusters/set_up_analysis.R*

 ## How do I run the code with the pre-prepared simulated data?
1. Enter the 'SI-from-clusters' folder and open the R project *.RProj* file. 

2. Open and run *Master_analysis.R*

 The transmission cloud for the simulated data has already been built and is included in the *data_simulated* folder. If you want to see how this was generated, take a look at *SI-from-clusters/Create_transcloud_sim.R* first.


 ## How can I run this with my own data?
 1. Create a new data folder containing 
 * a .fasta file of all viral sequences
 * a .csv file with 3 columns: sample_id = name of sequences in .fasta,	onset_date = case symptom onset date,	cluster_id = identifier of which cluster the case is in (can be character or numeric)

2. Enter the 'SI-from-clusters' folder and run *Create_transcloud_sim*, remembering to set your file paths and desired criteria for plausible transmission pairs. Proceed to *Master_analysis.R* to run the serial interval estimation. 

## Victorian SARS-CoV-2 data 
This repository also contains source data for the above manuscript, for cases of COVID-19 in Victoria, Australia. The *data_Victoria* folder contains GISAID accession numbers to access the sequences, and symptom onset times. The *SI-from-clusters-Victoria* folder then contains a version of the model code specifically adapted for the Victorian analysis.

## Authors
Jessica Stockdale (maintainer. jessica_stockdale@sfu.ca) <br />
Kurnia Susvitasari <br />
Anders Gonçalves da Silva <br />
Nicola Mulberry

 

