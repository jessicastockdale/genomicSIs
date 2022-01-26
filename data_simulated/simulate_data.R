###########################################################################
###
###  Simulate a outbreaks in several clusters with 'seedy' package
###
###########################################################################


library(seedy)
library(tidyverse)
library(seqinr)


################
# We need an initial genome to set off the simulation. We use the Wuhan reference genome of SARS-CoV-2, 
# but remove all non-CGAT bases
refgen <- ape::read.FASTA("../data_simulated/Wuhan_genome.fasta")
refgen <-as.numeric(refgen$MN908947.3[refgen$MN908947.3!="f0"])
# seedy needs sequences coded as 1,2,3,4
refgen <- refgen  %>% ifelse(. == 24, 4, .) %>% ifelse(. == 40, 1, .) %>% ifelse(. == 72, 3, .) %>% ifelse(. == 136, 2, .)
  

################
# Simulate the outbreaks
all_seqs <- vector(mode = "list", length = 0)
all_dat <- data.frame(sample_id = character(), onset_date = numeric(), cluster_id=numeric())
for (out in 1:2){
OB <- simulateoutbreak(init.sus=100, inf.rate = 0.22, rem.rate = 0.1, mut.rate = 0.1,
                 equi.pop = 2000, shape=flat, init.inf = 1, inoc.size = 1, 
                 samples.per.time = 1, samp.schedule = "random", 
                 full=FALSE, mincases = 30, ref.strain = refgen)



seqs <- librtoDNA(sampleID=OB$librstrains, libr=OB$libr, nuc=OB$nuc, 
               ref.strain=refgen, key=OB$librstrains, strings=TRUE)
names(seqs) <- OB$librstrains


#################
# Want to output: 1) .fasta with all sequences in
#                 2) .csv with metadata (sequence name, cluster, symptom onset)

# Fasta. Gather sequence for each infected case
total_seqs <- vector(mode = "list", length = nrow(OB$sampledata))
names(total_seqs) <- paste0("MYID_C",sprintf("%02d",out), OB$sampledata[,1])
for (i in 1:length(total_seqs)){
  total_seqs[[i]] <- seqs[names(seqs)==OB$sampledata[i,3]]
}

all_seqs <- c(all_seqs, total_seqs)

# CSV. Sequence metadata. Use infection times as symptom onset.
total_dat <- data.frame(sample_id = paste0("MYID_C",sprintf("%02d",out), OB$epidata[,1]), onset_date = OB$epidata[,2], cluster_id=c(rep(out, length(OB$epidata[,2]))))

all_dat <- rbind(all_dat, total_dat)

}
write.fasta(all_seqs, names(all_seqs), file.out = "../data_simulated/sim_genomes.fasta")
# Convert onset column to dates (make up origin date 01/04/20)
all_dat[,2] <- as.Date(all_dat[,2], origin="2020-04-01", format="%Y-%m-%d")
write.csv(all_dat, "../data_simulated/sim_metadata.csv")

