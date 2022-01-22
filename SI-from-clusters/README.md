# Estimating serial intervals from genomic clusters 

**Master_analysis.R**: main file from which to run serial interval estimation. This calls the below scripts

**set_up_analysis.R**: Called from master_analysis. Reads in data, libraries etc. 

**run_analysis.R**: Called from master_analysis. Performs the serial interval estimation

**results_compile.R**: Called from master_analysis. Generates plots and results table for serial interval estimation.

**Figs_epicurves.R**: Generates epi curves for each cluster, by onset and collection date. Called from master_analysis.
 
**Figs_clevelanddot.R**: Generates cleveland dot plots for each cluster, showing onset and collection dates of all cases. Called from master_analysis.

**Figs_trees.R**: Generates plots of the individual trees sampled in the serial interval analysis. Called from master_analysis.

**Figs_treenetworks.R**: Generates network plots of the space of trees sampled in the serial interval analysis. Called from master_analysis.

**estimate_Rt .R**: Generates Rt estimates from the cluster-specific serial interval estimates. Called from master_analysis.

**Create_transcloud.R**: Uses the transtreesampler package to build the transmission cloud. This needs to be run before Master_analysis.

**sensitivityanalysis_prior.R**: Performs sensitivity analysis to prior distributions for pi and w.

**pub_restable.R**: Export results table generated in Master_analysis to a csv (used to generate LaTeX tables in manuscript)


