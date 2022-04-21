###########################################################################
### In this script we tie together all the code for estimating serial 
### intervals in case clusters, for multiple waves of data.
###
### To run this, you will need:
### a transmission cloud: all possible infector-infectee pairs.
###    The transmission cloud must have the following column names:
###    cluster_id, case_i, case_j,  distance, onset_i, onset_j, onset_diff.
###    See Create_transcloud.R 
###
###  To adapt this analysis for a different dataset, you will need to edit the
###  set_up_analysis.R script with new filenames, cluster names etc.
###  The remaining scripts should then run as-is.
###
###########################################################################


###########################################################################
### Set up - select various options 
# Which wave of data do you want to work with?
which.wave <- 1  # set as wave 1 or wave 2
# Do you want to run the analysis with coprimary transmission in the model?
coprim.transm <- TRUE # Set TRUE if you want to use the model including coprimary, 
                      # or FALSE without (i.e. just compound geometric gamma)
# Do you want to include a prior on pi?
pi.model <- "prior" # Choices are "prior" and "none"
# Do you also want to include a prior on w?
w.model <- "prior" # Choices are "prior" and "none"
# Do you also want to pool all the trees from all the clusters together?
pool.trees <- TRUE # Set TRUE if you want to pool clusters for a single estimate as well as running separately
# Do you want to relabel the clusters when we output to figures (to match manuscript text)?
relabel=FALSE

# If you want to put a prior on pi, what beta prior parameters should we use?
if (pi.model=="prior"){pi.info <- c(beta_a = 12, beta_b = 11)
                        curve(dbeta(x,pi.info[1],pi.info[2]), ylab="density", main="Beta prior distribution for pi") }
# Otherwise, we don't need to decide anything about pi
if (pi.model=="none"){pi.info <- NULL}

# If you want to put a prior on w, what beta prior parameters should we use?
# NOTE: there is no w if coprim.transm=FALSE
if (coprim.transm == TRUE){
  if (w.model=="prior"){w.info <- c(beta_a = 12, beta_b = 11) 
  curve(dbeta(x,w.info[1],w.info[2]), ylab="density", main="Beta prior distribution for w") }
}
# Otherwise, we don't need to decide anything about w
if (w.model=="none"){w.info <- NULL}

# How many trees should we sample from each cluster? Must be at least 10
how.many = 10

###########################################################################
### Set-up - run the set-up script: loads libraries, data files, sets seed etc. 
source("set_up_analysis.R")


###########################################################################
### Run the serial interval analysis (also saves sampled trees and results to .rds file)
source("run_analysis.R")
res <- SI_estimation(data, names, coprim.transm, pi.model, w.model, pool.trees, pi.info, w.info, how.many, which.wave, progress.bars=TRUE)
res


# What is the chosen model?
model <- if_else(coprim.transm==FALSE, "nc", "cop")
if (pi.model=="prior"){model <- paste0(model, "_priorpi")}
if (w.model=="prior"){model <- paste0(model, "w")}
model <- paste0(model, "_w", which.wave)
# Save result object to a file
saveRDS(res, file = paste0(model, "_resultsobject.rda"))
# Read in a result object
res <- readRDS(file = paste0(model, "_resultsobject.rda"))





###########################################################################
### Read in metadata
metadata <- as_tibble(read.csv(paste0("../data_Victoria/epidata_w", which.wave, "clusters.csv"))) %>%
  # Rename columns to required format. Need: cluster_id, case_id, onset_date 
  rename(case_id = GISAID_ID) %>% 
  # Ensure onset_date is of date type
  mutate(onset_date=as.Date(onset_date, format = "%Y-%m-%d"))


# remove any cases with missing onset
metadata <- metadata[!is.na(metadata$onset_date),]


###########################################################################
### Plots of sampled tree networks

## sampled tree networks
source("Figs_treenetworks.R")
# What is the minimum posterior probability of a branch that should be shown in the networks?
# E.g. what is the minimum proportion of trees that this branch occurred in. 
threshold <- 0.1
plots_n <- treenets(metadata, names, coprim.transm, pi.model, w.model, threshold, which.wave)
plots_n
# The treenets function saves these networks as html. Should return TRUE for each cluster.


## Plot all of the individual sampled trees
source("Figs_trees.R")
plots_s <- samptrees(metadata, names, coprim.transm, pi.model, w.model, which.wave)
# These are saved these into their own subfolder within Figures




###########################################################################
### Relabel clusters, if required
which.names=1:length(names)
if (relabel){
  # Wave 1 = A, wave 2 = B etc. Numeric order by earliest onset date
  ord <- metadata %>% group_by(cluster_id) %>% slice(which.min(as.Date(onset_date))) %>% arrange(as.Date(onset_date)) %>% pull(cluster_id) # time order
  new.names <- paste0("cluster_", LETTERS[which.wave], 1:length(names))[match(sub(".*_", "", names), ord)]
  lookup <- tibble(cluster_id = sub(".*_", "", names), new.names = sub(".*_", "", new.names))
  data <- data %>% mutate(cluster_id = as.character(cluster_id)) %>% left_join(lookup, by = "cluster_id") %>% select(-cluster_id) %>% rename(cluster_id=new.names)
  metadata <- metadata %>% mutate(cluster_id = as.character(cluster_id)) %>% left_join(lookup, by = "cluster_id") %>% select(-cluster_id) %>% rename(cluster_id=new.names)
  names <- new.names
}else{lookup=NA}

  
###########################################################################
### Generate and print results table and figures (saves all to pdf)

source("results_compile.R")
results_gen <- SI_results(res, data, names, coprim.transm, pi.model, w.model, pool.trees, which.wave, which.names, pi.info, w.info)
grid.arrange(grobs=results_gen$Table)


if (which.wave==2){ # Figures by exposure site type (wave 2 only)
  source("Figs_exposure.R")
  plots_exposure <- exposure_types(res, metadata, names, coprim.transm, pi.model, w.model, pool.trees, which.wave, which.names)
  plots_exposure
}


###########################################################################
### Epi curves of each cluster
#
source("Figs_epicurves.R")
plots_e <- epi_curves(metadata, which.wave)
plots_e




###########################################################################
### Estimate Rt from the cluster-specific serial intervals
# NOTE: this requires original metadata to run
source("estimate_Rt.R")
rt <- estimate.rt(metadata, coprim.transm, pi.model, w.model, pool.trees, which.wave)

###########################################################################






