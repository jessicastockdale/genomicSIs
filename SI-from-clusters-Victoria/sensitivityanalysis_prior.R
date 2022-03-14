###########################################################################
### In this script we perform a sensitivity analysis to the prior
### distributions for pi and w
###
###########################################################################

library(tidyverse)
library(RColorBrewer)
Scenarios <- tibble(shape1 = c(17, 8, 3.5, 42, 12), shape2 = c(11,11,3,40,11))
Scen_names <- c("Baseline", "Increased mean", "Decreased mean", "Increased standard deviation", "Decreased standard deviation")
# 1. larger mean, 2. smaller mean, 3. larger sd, 4. smaller sd, 5. baseline

# Plot the scenarios
x <- seq(0,1,length.out=1000)
scurves <- tibble(x = x, Baseline = dbeta(x,12,11), "Increased mean" = dbeta(x,17,11), "Decreased mean" = dbeta(x,8,11),
                  "Increased standard deviation" = dbeta(x,3.5,3), "Decreased standard deviation" = dbeta(x,42,40))


mycols <- colorRampPalette(brewer.pal(8, "Spectral"))(4)
mycols <- c("#000000", mycols)

scurves %>% pivot_longer(!x, names_to = "Prior scenario") %>% 
  mutate(`Prior scenario` = factor(`Prior scenario`, levels=c("Baseline", "Increased mean","Decreased mean","Increased standard deviation","Decreased standard deviation"))) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_line(aes(color = `Prior scenario`), lwd=0.7) + theme_classic() + theme(legend.position="top", text = element_text(size=14)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  xlab(expression('Prior probability (' ~pi*' or '~w*')')) + ylab("Density") + scale_colour_manual(values = mycols)
ggsave("../Figures/prior_scenarios.pdf")



### Set up - select various options 
# Which wave of data do you want to work with?
which.wave <- 1  # set as wave 1 or wave 2
# Do you want to run the analysis with coprimary transmission in the model?
coprim.transm <- TRUE # Set TRUE if you want to use the model including coprimary
# Do you want to include a prior on pi, hold pi fixed or have unrestricted pi?
pi.model <- "prior" # Choices are "prior" and "none"
# Do you also want to include a prior on w, or have unrestricted w?
w.model <- "prior" # Choices are "prior" and "none"
# Do you also want to pool all the trees from all the clusters together?
pool.trees <- FALSE # Set TRUE if you want to pool clusters for a single estimate as well as running separately
# Do you want to relabel the clusters when output to figures?
relabel=FALSE
# How many trees should we sample from each cluster? Must be at least 10
how.many = 100


store_res <- vector(mode = "list", length = 5)
for (scen in 1:5){
  
  ###########################################################################
  # If you want to put a prior on pi, what beta prior parameters should we use?
  if (pi.model=="prior"){pi.info <- c(beta_a = Scenarios$shape1[scen], beta_b = Scenarios$shape2[scen])#(beta_a = 12, beta_b = 11)
  curve(dbeta(x,pi.info[1],pi.info[2]), ylab="density", main="Beta prior distribution for pi") }
  # Otherwise, we don't need to decide anything about pi
  if (pi.model=="none"){pi.info <- NULL}
  
  # If you want to put a prior on w, what beta prior parameters should we use?
  # NOTE: there is no w if coprim.transm=FALSE
  if (coprim.transm == TRUE){
    if (w.model=="prior"){w.info <- c(beta_a = Scenarios$shape1[scen], beta_b = Scenarios$shape2[scen]) # beta_a = 12, beta_b = 11)
    curve(dbeta(x,w.info[1],w.info[2]), ylab="density", main="Beta prior distribution for w") }
  }
  # Otherwise, we don't need to decide anything about w
  if (w.model=="none"){w.info <- NULL}
  
  
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
  if (which.wave==1){
    # wave 1 metadata
    metadata <- as_tibble(read.csv("../data_simulated/epidata_simclusters.csv")) %>%
      # Rename columns to required format. Need: cluster_id, case_id, onset_date (whether this be symptom, diagnosis)
      rename(case_id = sample_id) %>%
      # Ensure onset_date is of date type
      mutate(onset_date=as.Date(onset_date, format = "%Y-%m-%d"))
  }else{
    print("There is no wave 2 simulated data")
  }
  # remove any cases with missing symptom onset
  metadata <- metadata[!is.na(metadata$onset_date),]
  
  
  
  
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
    
  }
  
  ###########################################################################
  ### Generate and print results table and figures (also saves all to pdf)
  
  source("results_compile.R")
  results_gen <- SI_results(res, data, names, coprim.transm, pi.model, w.model, pool.trees, which.wave, which.names)
  grid.arrange(grobs=results_gen$Table)
  
  store_res[[scen]] <-res
  
}

# Save result object to a file
saveRDS(store_res, file = paste0(model, "_simulationresults.rda"))
# Read in a result object
store_res <- readRDS(file = paste0(model, "_simulationresults.rda"))




### Figures across scenarios

# Bars
names_formatted <- substring(names, 9)
res_dat <- tibble(cluster = rep(names_formatted, 5),
                  Scenario = rep(Scen_names[c(2:5,1)], each=length(names)),
                  mu = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$par[1]) })) })  ),
                  mu_ciL = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[1,1]) })) })  ),
                  mu_ciU = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[1,2]) })) })  ),
                  sigma = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$par[2]) })) })  ),
                  sigma_ciL = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[2,1]) })) })  ),
                  sigma_ciU = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[2,2]) })) })  ),
                  pi = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$par[3]) })) })  ),
                  pi_ciL = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[3,1]) })) })  ),
                  pi_ciU = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[3,2]) })) })  )
)
if (coprim.transm == TRUE){
  res_dat$w = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$par[4]) })) })  )
  res_dat$w_ciL = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[4,1]) })) })  )
  res_dat$w_ciU = unlist(lapply(store_res, function(olist){ unlist(lapply(olist$results[which.names], function(list) { return(list$estimates$confint[4,2]) })) })  )
}

mycols <- colorRampPalette(brewer.pal(8, "Spectral"))(4)
mycols <- c("#000000", mycols)
res_dat$cluster <- factor(res_dat$cluster, levels = c(mixedsort(names_formatted), "ALL"))
res_dat$Scenario <- factor(res_dat$Scenario, levels = Scen_names)


p1 <- ggplot(res_dat, aes(x=cluster, y=mu, col=Scenario, group=Scenario)) +
  geom_errorbar(aes(ymin = mu_ciL, ymax = mu_ciU), width=0.2, position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5))+ ylim(0,max(11, max(res_dat$mu_ciU))) +
  labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu)) 
p1 <- p1 + theme_classic() + scale_colour_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14),
        legend.position="top", legend.direction="horizontal") + guides(colour=guide_legend(nrow=2,byrow=TRUE))

legend <- cowplot::get_legend(p1)
grid.newpage()
pdf(file = paste0("../Figures/sensitivity_legend.pdf"),  width= 8, height=1.5)
grid.draw(legend)
dev.off()


p1 <- ggplot(res_dat, aes(x=cluster, y=mu, col=Scenario, group=Scenario)) +
  geom_errorbar(aes(ymin = mu_ciL, ymax = mu_ciU), width=0.2, position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5))+ ylim(0,max(11, max(res_dat$mu_ciU))) +
  labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu)) 
p1 <- p1 + theme_classic() + scale_colour_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14),
        legend.position="none")
pdf(file = paste0("../Figures/sensitivity_mubars.pdf"),  width= 8, height=2.8)
  p1
dev.off()




p2 <- ggplot(res_dat, aes(x=cluster, y=sigma, col=Scenario, group=Scenario)) +
  geom_errorbar(aes(ymin = sigma_ciL, ymax = sigma_ciU), width=0.2, position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5))+ ylim(0,max(7, max(res_dat$sigma_ciU))) +
  labs(title=paste0(measure, " Standard Deviation"),x="Cluster", y = expression(sigma)) 
p2 <- p2 + theme_classic() + scale_colour_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14),
        legend.position="none") 
pdf(file = paste0("../Figures/sensitivity_sigmabars.pdf"),  width= 8, height=2.8)
p2
dev.off()


p3 <- ggplot(res_dat, aes(x=cluster, y=pi, col=Scenario, group=Scenario)) +
  geom_errorbar(aes(ymin = pi_ciL, ymax = pi_ciU), width=0.2, position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5)) + ylim(0,1) +
  labs(title="Sampling Proportion",x="Cluster", y = expression(pi)) 
p3 <- p3 + theme_classic() + scale_colour_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14),
        legend.position="none") 
pdf(file = paste0("../Figures/sensitivity_pibars.pdf"),  width= 8, height=2.8)
p3
dev.off()


p4 <- ggplot(res_dat, aes(x=cluster, y=w, col=Scenario, group=Scenario)) +
  geom_errorbar(aes(ymin = w_ciL, ymax = w_ciU), width=0.2, position=position_dodge(0.5))+
  geom_point(position=position_dodge(0.5)) + ylim(0,1) +
  labs(title="Proportion Non-coprimary",x="Cluster", y = "w") 
p4 <- p4 + theme_classic() + scale_colour_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14),
        legend.position="none") 
pdf(file = paste0("../Figures/sensitivity_wbars.pdf"),  width= 8, height=2.8)
p4
dev.off()


