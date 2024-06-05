#####################################################################
### This script creates the 'transmission cloud', for simulated data
###
### That is, a list of all plausible infector/infectee pairs who meet 
### a set of genomic and symptom onset time criteria
###
### It reads in the alignment (.fasta) and metadata, and
### outputs the transmission cloud as .csv file in which each row is a
### potential infector-infectee pair.
###
#####################################################################

# Only need to run this the first time to install 'transtreesampler' package:
install.packages("remotes")
remotes::install_github("jessicastockdale/transtreesampler")
#

library(transtreesampler)
library(magrittr)

## prepare data
data_dir <- "../data_simulated"
fasta <- "sim_genomes.fasta"
metadata <- "sim_metadata.csv"

# load data ---------------------------------------------------------------
aln <- ape::read.FASTA(file.path(data_dir, fasta))
names(aln) <- sub('\\|.*', '', labels(aln))
mdt <- readr::read_csv(file.path(data_dir, metadata))


# set params --------------------------------------------------------------
# Do you want to get serial intervals or collection/diagnosis intervals?
which.type <- quote(onset_date) # or quote(diagnosis_date) 

cluster_col <- quote(cluster_id) # which column indicates cluster names
sample_id_col <- quote(sample_id) # which column indicates sample ids
onset_date_col <- which.type
min_cluster_size <- 15
max_genetic_dist <- 1.1/29306 # approx. 1 SNP equivalent
min_day_diff <- 1
max_day_diff <- 35
output_filename <- "../data_simulated/sim_transmission_cloud.csv"

# remove cases with no symptom onset
mdt <- mdt[!is.na(mdt$onset_date),]  


# transmission cloud ------------------------------------------------------

(tmp <- transtreesampler::id_transmission_pairs(mdt,
                                                aln,
                                                max_genetic_dist = max_genetic_dist,
                                                min_day_diff = min_day_diff,
                                                max_day_diff = max_day_diff,
                                                onset_date_col = onset_date_col,
                                                case_id_col = sample_id_col,
                                                cluster_col = cluster_col, min_cluster_size = min_cluster_size))

tmp %>%
  dplyr::select(cluster_id, genetic_dists) %>%
  tidyr::unnest(genetic_dists) %>%
  readr::write_csv(output_filename)

# Create csv with relevant clusters only
mdt <- mdt[mdt$cluster_id %in% tmp$cluster_id,]  
write.csv(mdt, "../data_simulated/epidata_simclusters.csv")

#print list of case IDs in clusters
dir.create(file.path("../data_simulated/cluster_ids"), showWarnings = FALSE)
for(i in 1:length(tmp$cluster_id)){ 
  x <- unique(tmp$data[[i]]$sample_id)
  filename <- file.path("../data_simulated/cluster_ids",paste("cluster_", tmp$cluster_id[[i]],"_IDs.txt",sep=""))
  write.table(x, file = filename, row.names=FALSE, col.names=FALSE)
}
