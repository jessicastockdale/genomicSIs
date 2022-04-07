## This code sets up the analysis:
# sets the seed,
# loads all required packages,
# reads in the data

## Set seed
set.seed(58392)

## Load libraries

# Package names
packages <- c("ggplot2", "tidyverse", "gridExtra", "gridBase", "reshape2", "lattice",
              "future", "RColorBrewer", "stringr", "igraph", "visNetwork", "numDeriv",
              "doParallel", "doSNOW", "tcltk", "EpiEstim", "ggpubr", "cowplot", "scales",
              "padr", "gtools", "grid", "ape", "seqinr", "HDInterval")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
rm(installed_packages, packages)

# Create directory for sampled trees
dir.create(file.path("sampled_trees"), showWarnings = FALSE)

## Read in the data
data <- as_tibble(read.csv(paste0("../data_Victoria/w", which.wave, "_transmission_cloud.csv")))
# Ensure cluster column is called cluster_id
data <- rename(data, cluster_id = cluster_id)
# Set cluster names
names <- sort(paste0("cluster_", unique(data$cluster_id)))
for(i in names){ 
  filepath <- file.path("../data_Victoria/cluster_ids",paste(i,"_IDs.txt",sep=""))
  assign(i, read.table(filepath))
}


