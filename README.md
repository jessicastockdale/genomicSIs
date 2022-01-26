# Estimation of serial intervals using viral sequences

This repository contains code for estimating the serial interval distribution in case clusters using genomic data, taking incomplete sampling into account.

This is supporting code for <br />
**'Genomic epidemiology offers high resolution estimates of serial intervals for COVID-19'** <br />
Jessica E. Stockdale, Kurnia Susvitasari, Paul Tupper, Benjamin Sobkowiak, Nicola Mulberry, Anders Gonçalves da Silva, Anne E. Watt, Norelle Sherry, Corinna Minko, Benjamin P. Howden, Courtney R. Lane, and Caroline Colijn.

## File/folder contents:
* 'SI-from-clusters' contains all user-facing code and the R project file. A readme with full details of all scripts is included in this folder.
* 'SI with coprimary' contains utility functions for the underlying statistical model, including coprimary transmission.
* 'SI with noncoprimary' contains equivalent utilities for a model with no coprimary transmission.
* 'data_simulated' contains a simulated dataset which is ready to be used with the above code
 
 All code was prepared and tested under R version 4.1.0 (2021-05-18). Package dependencies are listed and will be automatically installed in * *SI-from-clusters/set_up_analysis.R* *
 

