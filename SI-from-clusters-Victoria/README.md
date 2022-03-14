# Estimating serial intervals in clusters of cases using viral sequences and symptom onset times, taking incomplete sampling into account.

This folder is set up to run the analysis for two waves of Victorian COVID-19 data, as in the 'data_Victoria' folder

**genomicSIs_Victoria.RProj**: R project file.

**Create_transcloud_w1.R & Create_transcloud_w2.R**: Uses the [transtreesampler](https://github.com/andersgs/transtreesampler) package to build the transmission cloud. This is the first script you should run.


**Master_analysis.R**: main file from which to run serial interval estimation. This calls the below scripts:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **set_up_analysis.R**: Called from master_analysis. Reads in data, libraries etc. 

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **run_analysis.R**: Called from master_analysis. Performs the serial interval estimation

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **results_compile.R**: Called from master_analysis. Generates plots and results table.

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Figs_epicurves.R**: Generates epi curves for each cluster, by onset date. Called from master_analysis.

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Figs_trees.R**: Generates plots of the individual trees sampled in the serial interval analysis. Called from master_analysis.

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **Figs_treenetworks.R**: Generates network plots of the space of trees sampled in the serial interval analysis. Called from master_analysis.

 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; **estimate_Rt .R**: Generates Rt estimates from the cluster-specific serial interval estimates. Called from master_analysis.


**sensitivityanalysis_prior.R**: Performs sensitivity analysis to prior distributions for pi and w.

**pub_restable.R**: Export results table generated in Master_analysis to a csv (used to generate LaTeX tables in manuscript)


