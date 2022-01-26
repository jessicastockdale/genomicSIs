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
              "padr", "gtools", "grid", "ape", "seqinr")
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
if (which.wave==1){data <- as_tibble(read.csv("../data_simulated/sim_transmission_cloud.csv"))
# Ensure cluster column is called cluster_id
data <- rename(data, cluster_id = cluster_id)
# Set cluster names
names <- sort(paste0("cluster_", unique(data$cluster_id)))
for(i in names){ 
  filepath <- file.path("../data_simulated/cluster_ids",paste(i,"_IDs.txt",sep=""))
  assign(i, read.table(filepath))
}
}else{ # wave 2:
  print("Simulated dataset only has 1 wave")
}


