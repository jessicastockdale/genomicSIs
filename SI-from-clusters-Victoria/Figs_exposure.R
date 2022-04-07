
########################################################
### Create figures categorized by exposure site type ###
########################################################


exposure_types <- function(results, metadata, names, coprim.transm=TRUE, pi.model="none", w.model="none", pool.trees=FALSE, which.wave = 1, which.names = 1:length(names)){
  

  ## Set up model name, for file naming system
  if (coprim.transm==FALSE){model <- "nc_"
  }else{                    model <- "cop_"
  }
  if (pi.model=="prior"){
    model <- paste0(model, "priorpi")
  }
  if (w.model=="prior"){
    model <- paste0(model, "w")
  }
  if (pool.trees==TRUE){
    model <- paste0("pooled_", model)
  }
  model <- paste0(model, "_w", which.wave)
  measure <- "Serial Interval"
  print(paste0("Confirming, chosen model is: ", model)) 
  # Create directory for model results, if it doesn't already exist
  dir.create(file.path(paste0("../Figures/", model)), showWarnings = FALSE)
  

## Plot the mean and CI for each cluster - using errorbars, coloured by exposure site type
names_formatted <- substring(names, 9)
# Of pooled trees
if (pool.trees==TRUE){
  res_dat_pooled <- tibble(cluster = rep("All", 1),
                    mu = results$results_pooled$estimates$par[1],
                    mu_ciL = results$results_pooled$estimates$confint[1,1],
                    mu_ciU = results$results_pooled$estimates$confint[1,2],
                    sigma = results$results_pooled$estimates$par[2],
                    sigma_ciL = results$results_pooled$estimates$confint[2,1],
                    sigma_ciU = results$results_pooled$estimates$confint[2,2],
                    pi = results$results_pooled$estimates$par[3],
                    pi_ciL = results$results_pooled$estimates$confint[3,1],
                    pi_ciU =results$results_pooled$estimates$confint[3,2],
  )
  if (coprim.transm == TRUE){
    res_dat_pooled$w = results$results_pooled$estimates$par[4]
    res_dat_pooled$w_ciL = results$results_pooled$estimates$confint[4,1]
    res_dat_pooled$w_ciU = results$results_pooled$estimates$confint[4,2]
  }
}

# Now in each cluster
res_dat <- tibble(cluster = names_formatted,
                  mu = unlist(lapply(results$results, function(list) { return(list$estimates$par[1]) })),
                  mu_ciL = unlist(lapply(results$results, function(list) { return(list$estimates$confint[1,1])})),
                  mu_ciU = unlist(lapply(results$results, function(list) { return(list$estimates$confint[1,2])})),
                  sigma = unlist(lapply(results$results, function(list) { return(list$estimates$par[2]) })),
                  sigma_ciL = unlist(lapply(results$results, function(list) { return(list$estimates$confint[2,1])})),
                  sigma_ciU = unlist(lapply(results$results, function(list) { return(list$estimates$confint[2,2])})),
                  pi = unlist(lapply(results$results, function(list) { return(list$estimates$par[3]) })),
                  pi_ciL = unlist(lapply(results$results, function(list) { return(list$estimates$confint[3,1])})),
                  pi_ciU = unlist(lapply(results$results, function(list) { return(list$estimates$confint[3,2])})),
)
if (coprim.transm == TRUE){
  res_dat$w = unlist(lapply(results$results, function(list) { return(list$estimates$par[4]) }))
  res_dat$w_ciL = unlist(lapply(results$results, function(list) { return(list$estimates$confint[4,1])}))
  res_dat$w_ciU = unlist(lapply(results$results, function(list) { return(list$estimates$confint[4,2])}))
}

# Add exposure site types
res_dat <- metadata[,c("cluster_id", "exposure_site_category")] %>% distinct() %>%  
                  merge(res_dat, by.y = c("cluster"), by.x = c("cluster_id"), all.y=TRUE ) %>% rename("cluster"="cluster_id")


# Add pooled trees if available
ncols <- length(unique(res_dat$exposure_site_category)) 
if (pool.trees==TRUE){
  ncols <- ncols + 1
  res_dat_pooled$exposure_site_category <- "All"
  res_dat <- rbind(res_dat, res_dat_pooled)
}
mycols <- brewer.pal(8, "Set1")[-6]
res_dat$cluster <- factor(res_dat$cluster, levels = c(sort(names_formatted), "All"))

# sort clusters by ...
levs <- res_dat$cluster[order(res_dat$mu)]
levs <- c(levs[1:(which(levs=="All")-1)], levs[(which(levs=="All")+1):length(levs)], levs[which(levs=="All")])
res_dat$cluster <- factor(res_dat$cluster, levels=levs)

p1 <- ggplot(res_dat, aes(x=cluster, y=mu, col=exposure_site_category)) +
  geom_errorbar(aes(ymin = mu_ciL, ymax = mu_ciU), width=0.2)+
  geom_point()+
  labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu)) 
p1 <- p1  + theme_classic() + scale_colour_manual(values = mycols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14), legend.position=c(0.3, 0.85), legend.direction = "horizontal") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
p1

p2 <- ggplot(res_dat, aes(x=cluster, y=sigma, col=exposure_site_category)) +
  geom_errorbar(aes(ymin = sigma_ciL, ymax = sigma_ciU), width=0.2)+
  geom_point()+
  labs(title=paste0(measure, " SD"),x="Cluster", y = expression(sigma))
p2 <- p2  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

p3 <- ggplot(res_dat, aes(x=cluster, y=pi, col=exposure_site_category)) +
  geom_errorbar(aes(ymin = pi_ciL, ymax = pi_ciU), width=0.2)+
  geom_point()+ylim(0,1)+
  labs(title="Sampling Proportion",x="Cluster", y = expression(pi))
p3 <- p3  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

if (coprim.transm == TRUE){
  p6 <- ggplot(res_dat, aes(x=cluster, y=w, col=exposure_site_category)) + 
    geom_errorbar(aes(ymin = w_ciL, ymax = w_ciU), width=0.2)+
    geom_point()+ylim(0,1)+
    labs(title="Proportion non-coprimary",x="Cluster", y = "w")
  p6 <- p6  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          text = element_text(size=14)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}

# Print to pdf
pdf(file = paste0("../Figures/", model, "/bars_exptype.pdf"), height = 10, width= 11.7)
if (coprim.transm == FALSE){
  grid.arrange(p1, p2, p3, nrow=2)
}else{grid.arrange(p1, p2, p3, p6, nrow=4)}
dev.off()
if (coprim.transm == FALSE){
  grid.arrange(p1, p2, p3, nrow=2)
}else{grid.arrange(p1, p2, p3, p6, nrow=4)}





### Scatterplots of date vs mean SI and cluster size vs pi 

# remove 'all' row of res_dat, if there is one
if (pool.trees==TRUE){
  res_dat <- res_dat[-(nrow(res_dat)),]
}

# data structure
min_onset <- metadata %>% group_by(cluster_id) %>% summarise(minons = min(onset_date, na.rm=TRUE))
scatterdat <- tibble(cluster = res_dat$cluster, mu = res_dat$mu, pi = res_dat$pi, initial.onset = as.Date(min_onset$minons), 
                     cl.size = table(metadata$cluster_id), Type = res_dat$exposure_site_category)
# Change NA to "NA"
scatterdat$Type[is.na(scatterdat$Type)] <- "NA"

ncols <- length(unique(scatterdat$Type))
cols <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)

sp1 <- ggscatter(scatterdat, x = "initial.onset", y = "mu",
                 color = "Type", palette = cols) +
  labs(x="Initial symptom onset", y = expression(mu)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), text = element_text(size=16)
  ) + scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))

sp2 <- ggscatter(scatterdat, x = "initial.onset", y = "pi",
                 color = "Type", palette = cols) +
  labs(x="Initial symptom onset", y = expression(pi)) +
  theme(text = element_text(size=16)
  )  + scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))

sp3 <- ggscatter(scatterdat, x = "cl.size", y = "mu",
                 color = "Type", palette = cols) +
  labs(x="Cluster size", y = expression(mu))   +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(), text = element_text(size=16)
  )  + scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))   

sp4 <- ggscatter(scatterdat, x = "cl.size", y = "pi",
                 color = "Type", palette = cols) +
  labs(x="Cluster size", y = expression(pi)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(), text = element_text(size=16)
  )    + scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1))


# height = 8.3, width= 11.7
ggarrange(sp1, sp3, sp2, sp4, nrow=2, ncol=2, common.legend = TRUE, legend="top")

# Print to pdf
pdf(file = paste0("../Figures/", model, "/corrs_exptype.pdf"), height = 10, width= 11.7)
ggarrange(sp1, sp3, sp2, sp4, nrow=2, ncol=2, common.legend = TRUE, legend="top")
dev.off()



# Correlations
print(cor.test(as.numeric(scatterdat$initial.onset),scatterdat$mu))
print(cor.test(as.numeric(scatterdat$initial.onset),scatterdat$mu, method = "spearman"))

print(cor.test(as.numeric(scatterdat$initial.onset),scatterdat$pi))
print(cor.test(as.numeric(scatterdat$initial.onset),scatterdat$pi, method = "spearman"))

print(cor.test(as.numeric(scatterdat$cl.size),scatterdat$mu))
print(cor.test(as.numeric(scatterdat$cl.size),scatterdat$mu, method = "spearman"))

print(cor.test(as.numeric(scatterdat$cl.size),scatterdat$pi))
print(cor.test(as.numeric(scatterdat$cl.size),scatterdat$pi, method = "spearman"))

}

