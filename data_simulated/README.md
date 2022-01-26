# Simulated cluster outbreaks

This folder contains simulated data on 2 cluster outbreaks. These cluster outbreaks were simulated using the R package [seedy](https://doi.org/10.1371/journal.pone.0129745). These are formatted and ready to be input to the serial interval estimation procedure in the *SI-from-clusters* folder.

Cluster 1 - 43 cases
Cluster 2 - 46 cases

Both clusters simulated with: 
* population size 100 
* infection rate 0.22
* removal rate 0.1
* mutation rate 0.1
* sampling of each infected case at a random point during their infectious period.

## File contents:
* 'sim_genomes' contains the simulated viral sequence for each case
* 'sim_metadata.csv' contains the metadata of the simulated outbreaks. For each case, identified with a *sample_id*, this lists the symptom onset time *onset_date* and cluster membership *cluster_id*.
* 'epidata_simclusters.csv' contains the metadata of only those clusters selected for serial interval analysis e.g. with at least 15 cases. For this simulated data, that is all clusters and so epidata_simclusters.csv = sim_metadata.csv
* 'simulate_data.R' contains the code used to simulate the clusters with *seedy*
* 'Wuhan_genome.fasta' contains the SARS-CoV-2 Wuhan reference genome. This was used to intialize the sequences in the simulation.
* 'sim_transmission_cloud.csv' contains the transmission cloud for the simulated data, as output by *SI-from-clusters/Create_transcloud_sim.R*.
* 'cluster_ids' contains .txt files listing which cases are in which clusters. This is output by *SI-from-clusters/Create_transcloud_sim.R*.

 

 

