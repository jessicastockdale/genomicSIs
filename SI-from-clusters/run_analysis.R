## This code run the serial interval estimation, given the options chosen in
# Master_analysis.R

SI_estimation <- function(data, names, coprim.transm=TRUE, pi.model="none", w.model="none", pool.trees=FALSE, pi.info=NA, w.info=NA, how.many=10, which.wave=1, progress.bars=TRUE){
  
  # Read in the relevant model for the chosen model options  
  if (coprim.transm==FALSE){
    source("../SI with noncoprimary/multitreeSI.R")
    source("../SI with noncoprimary/cgg_utilities.R")
    model <- "nc"
  }else{
    source("../SI with coprimary/SIestim.R")
    source("../SI with coprimary/coputilities.R")
    model <- "cop"
  }
  
  if (pi.model=="prior"){
    model <- paste0(model, "_priorpi")
  }
  if (w.model=="prior"){
    model <- paste0(model, "w")
  }
  print(paste0("Confirming, chosen model is: ", model))  
  # Create directory for model results, if it doesn't already exist
  dir.create(file.path(paste0("../Figures/", model)), showWarnings = FALSE, recursive = TRUE)
  
  
  ## Configuration of the SI estimation
  init_ <- list(mu=8,sigma=2,pi=0.8) #initial values for alpha, beta, pi
  low_ <- c(0.001, 0.001, 0.3) #lower bound of each estimate
  upp_ <- c(30, 30, 1.001) #upper bound of each estimate
  cores_ <- availableCores() - 1 #how many core processors
  # initial, upper and lower values for w where needed
  if (coprim.transm==TRUE){init_$w <- 0.8
  low_ <- c(low_, 0.01)
  upp_ <- c(upp_, 1.0)} 
  
  # Collect all the cluster results:
  howmany.run <- round(how.many*1.2) # We run extra, to account for the fact that some may not converge
  
  
  ## Regular (non-pooled) tree version
  # run per-cluster
  results_sam1 <- vector(mode = "list", length = length(names))
  for (i in 1:length(results_sam1)){
    while (nrow(results_sam1[[i]]$record) != how.many || is.null(results_sam1[[i]]$record)){
      print(paste0("Currently running ", names[i]))
      
      trees_sam <- list()
      # Select only pairs in this cluster
      cluster_data <- data[data$cluster_id==sub(".*_", "", names[i]),]
      
      
      # Tree sampling
      # Uniformly at random among all possible infectors:
      #for (j in 1:howmany.run){
      # For each infectee, we select an infector at random, in each tree
      #  cluster_datar <- cluster_data[sample(1:nrow(cluster_data)),]
      #  trees_sam[[j]] <- cluster_datar[!duplicated(cluster_datar$case_j),]
      #}
      
      # Probabilistic sampling - sample according to:
      #                         (i) shorter genomic distance is more likely
      #                         (ii) shorter time is more likely (either
      #                               totally, or above 1 day).
      # Genetic and time components are normalized to have equal impact.
      for (j in 1:howmany.run){
        maxdist <- max(cluster_data$distance)
        maxtime <- max(cluster_data$onset_diff)
        if (maxdist==0 && maxtime==0){stop("All cases in the cluster have 0 genomic distance and 0 onset difference")}
        # option a: 1 day SI is less likely
        # option b: a more straightforward less-time-is-better
        trees_sam[[j]] <-  cluster_data[order(cluster_data$case_j),] %>% 
          group_by(case_j) %>%
          slice_sample(n=1, weight_by = 
                         # option a: 
                         # (if(maxtime==0){0}else{dgamma(onset_diff, shape=2, scale=5)/max(dgamma(onset_diff, shape=2, scale=5))}) + (if(maxdist==0){0}else{abs(distance-maxdist)/maxdist}) )
                         # option b:
                         pmax(0.00001, (if(maxtime==0){0}else{abs(onset_diff-maxtime)/maxtime}) + (if(maxdist==0){0}else{abs(distance-maxdist)/maxdist}) ))
        # If the probability is 0 or very small, set to 0.00001 to avoid errors
      }
      
      
      
      # Save sampled trees to file
      saveRDS(trees_sam, file = paste0("sampled_trees/",names[i], model, "_w", which.wave, ".rds"))
      
      # Plot SIs of 10 of the sampled trees
      pdf(paste0("../Figures/w", which.wave, "_obs_si_", names[i],".pdf"), height = 8.3, width= 11.7, paper="a4r")
      par(mfrow = c(2,5), oma = c(0, 0, 2, 0))
      for (j in 1:10){
        hist(trees_sam[[j]]$onset_diff,
             main = paste("Tree", j),
             xlab = "Observed serial interval (days)", col = "cadetblue")
      }
      mtext(paste("Observed Serial Intervals of 10 sampled trees in ", str_to_title(gsub("_", " ", names[i], fixed=TRUE))), outer = TRUE, cex = 1.2)
      dev.off()
      gc()
      
      
      #estimate the parameters using multitreeSI
      if (coprim.transm==FALSE){
        # Non-coprimary
        rec.unknownM <- multitreeSI(trees_sam, si.config=list(init=init_, lower=low_, upper=upp_, ci=0.95), nameSI="onset_diff", cores=cores_, howmany = how.many, pi_prior=pi.info, progress.bars=progress.bars)
      }else{
        # Coprimary
        rec.unknownM <- SIestim(trees_sam, si.config=list(par0=init_, lower=low_, upper=upp_, ci=0.95), method="OPTIM", nameSI="onset_diff", cores=cores_, howmany = how.many, pi_prior=pi.info, w_prior = w.info, progress.bars=progress.bars)
      }
      print(rec.unknownM$estimates)
      
      results_sam1[[i]] <- rec.unknownM
      rm(rec.unknownM)
      gc()
      if(nrow(results_sam1[[i]]$record) != how.many){
        print("Warning: too few converged, re-running with 50% more trees.")
        howmany.run <- round(howmany.run*1.5) # run more}
      }else{howmany.run <- round(how.many*1.2)}
    }
  }
  # Save results object to file
  saveRDS(results_sam1, file = paste0(model, "_w", which.wave, "_results.rds"))
  
  
  
  ## Pooled tree version
  if (pool.trees==TRUE){
    
    print("Running pooled clusters analysis")
    trees_sam <- vector(mode = "list", length = howmany.run*length(names))
    results_sam_p <- vector(mode = "list", length = 1)
    
    for (i in 1:length(names)){
      # Read in the first howmany.run sampled trees from the first part of the analysis, 
      # and then bind together in a single list
      trees_sam[((i-1)*howmany.run +1):(i*howmany.run)] <- readRDS(file = paste0("sampled_trees/",names[i], model, "_w", which.wave, ".rds"))[1:howmany.run]
    }
    
    #estimate the parameters 
    if (coprim.transm==FALSE){
      # Non-coprimary
      rec.unknownM <- multitreeSI(trees_sam, si.config=list(init=init_, lower=low_, upper=upp_, ci=0.95), nameSI="onset_diff", cores=cores_, howmany = how.many*length(names), pi_prior=pi.info, progress.bars=progress.bars)
    }else{
      # Coprimary
      rec.unknownM <- SIestim(trees_sam, si.config=list(par0=init_, lower=low_, upper=upp_, ci=0.95), method="OPTIM",nameSI="onset_diff", cores=cores_, howmany = how.many*length(names), pi_prior=pi.info, w_prior = w.info, progress.bars=progress.bars)
    }
    print(rec.unknownM$estimates)
    
    results_sam_p <- rec.unknownM
    
    # Save results object to file
    model <- paste0("pooled_", model)
    saveRDS(results_sam_p, file = paste0(model, "_w", which.wave, "_results.rds"))
    
    
    rm(rec.unknownM)
    gc()
    if(nrow(results_sam_p$record) != how.many*length(names)){
      print("Warning: too few pooled trees converged. Estimate will be taken from fewer trees than intended.")}
    
  }else{results_sam_p <- NA}
  
  
  
  return(list(results = results_sam1, results_pooled = results_sam_p))
}