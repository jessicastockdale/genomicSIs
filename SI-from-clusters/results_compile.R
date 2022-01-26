## This code creates various tables and figures from the SI estimation analysis, link to 
# this from Master_analysis.R

SI_results <- function(results, data, names, coprim.transm=TRUE, pi.model="none", w.model="none", pool.trees=FALSE, which.wave = 1, which.names = 1:length(names)){
  
 ## Histograms of observed serial intervals
  if (which.wave==1){
    pdf("../Figures/obs_si.pdf", height = 8.3, width= 11.7, paper="a4r")
  }else{pdf("../Figures/obs_si_w2.pdf", height = 8.3, width= 11.7, paper="a4r")} 
   
  par(mfrow = c(2,5), oma = c(0, 0, 2, 0))
  for (i in mixedsort(names)){
    hist(data$onset_diff[data$cluster_id==sub(".*_", "", i)],
         main = paste(str_to_title(gsub("_", " ", i, fixed=TRUE))),
         xlab = "Serial interval (days)", col = "cadetblue")
  }
  mtext("Transmission cloud: all observed serial intervals", outer = TRUE, cex = 1.2)
  dev.off()
  par(mfrow = c(2,5), oma = c(0, 0, 2, 0))
  for (i in mixedsort(names)){
    hist(data$onset_diff[data$cluster_id==sub(".*_", "", i)],
         main = paste(str_to_title(gsub("_", " ", i, fixed=TRUE))),
         xlab = "Serial interval (days)", col = "cadetblue")
  }
  mtext("Transmission cloud: all observed serial intervals", outer = TRUE, cex = 1.2)
  
  
  
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
  if (which.wave==2){
    model <- paste0(model, "_w2")
  }
  measure <- "Serial Interval"
  print(paste0("Confirming, chosen model is: ", model)) 
  # Create directory for model results, if it doesn't already exist
  dir.create(file.path(paste0("../Figures/", model)), showWarnings = FALSE)
  

 ## Results table
  results_table1 <- data.frame(results$results[[which.names[1]]]$estimates$par)
  for (i in 2:length(which.names)){
    results_table1[,i] <- results$results[[which.names[i]]]$estimates$par
  }
  names(results_table1) <- sub(".*_", "", names)
  results_table1 <- results_table1[ , mixedsort(names(results_table1))] # order the clusters
  # Add pooled trees column if available
  if (pool.trees==TRUE){
    results_table1 <- cbind(results_table1, data.frame(results$results_pooled$estimates$par))
    names(results_table1)[length(names)+1] <- c("ALL")
    }
  results_table1 <- results_table1 %>% 
    mutate_if(is.numeric, round, digits=2)
  if (coprim.transm==FALSE){
  row.names(results_table1) <- c("mu", "sigma", "pi")
  } else{row.names(results_table1) <- c("mu", "sigma", "pi","w")}
  
  # Save table to RDS
  saveRDS(results_table1, file = paste0("resultstable_", model, ".rds"))
  # Transform table into one of multiple grobs
  Table <- vector(mode = "list", length = ceiling(length(names)/5))
  for (i in 1:ceiling(length(names)/5)){
    Table[[i]] <- tableGrob(results_table1[,((i-1)*5+1):min(i*5, length(names))])
  }
  if (pool.trees==TRUE){Table[[i]] <- tableGrob(results_table1[,((i-1)*5+1):(min(i*5, length(names))+1)])}
  # Create the PDF 
  m1 <- marrangeGrob(grobs=Table,nrow=5,ncol=1)
  ggsave(paste0("../Figures/", model, "/results_table.pdf"), m1, height = 8.3, width= 11.7, paper="a4r")
  
  
 ## Histograms of cluster estimates
  nr <- ifelse(coprim.transm==FALSE, 3, 4)
  pdf(file = paste0("../Figures/", model, "/hists.pdf"), height = 8.3, width= 11.7, paper="a4r")
  par(mfcol = c(nr,5), oma = c(0, 0, 2, 0))
  for (i in mixedorder(names)){
    hist(results$results[[i]]$record$mu, xlab="mu", breaks=15, main=paste0("Cluster ", sub(".*_", "", names[i])))
    hist(results$results[[i]]$record$sigma, xlab="sigma", breaks=15, main=paste0("Cluster ", sub(".*_", "", names[i])))
    if (pi.model!="fixed"){hist(results$results[[i]]$record$pi, xlab="pi", breaks=15, main=paste0("Cluster ", sub(".*_", "", names[i])))
      }else{hist(results$results[[i]]$estimates$par[3], xlab="pi", breaks=15, main=paste0("Cluster ", sub(".*_", "", names[i])))}
   if(coprim.transm==TRUE){
     hist(results$results[[i]]$record$w, xlab="w", breaks=15, main=paste0("Cluster ", sub(".*_", "", names[i])))
   }
  }
  dev.off()
  
  
 ## Histograms of pooled estimates
  if (pool.trees==TRUE){
  if (coprim.transm==FALSE){
    pdf(file = paste0("../Figures/", model, "/hists_pooled.pdf"), height = 8.3, width= 11.7, paper="a4r")
    layout(matrix(1:3, nrow=3), widths=rep.int(1,3), heights=rep.int(1,1), respect=F)
    hist(results$results_pooled$record$mu, xlab="mu", breaks=15, main="All Clusters")
    hist(results$results_pooled$record$sigma, xlab="sigma", breaks=15, main="All Clusters")
    hist(results$results_pooled$record$pi, xlab="pi", breaks=15, main="All Clusters")
    dev.off()
  } else{pdf(file = paste0("../Figures/", model, "/hists_pooled.pdf"), height = 8.3, width= 11.7, paper="a4r")
    layout(matrix(1:4, nrow=4), widths=rep.int(1,3), heights=rep.int(1,1), respect=F)
    hist(results$results_pooled$record$mu, xlab="mu", breaks=15, main="All Clusters")
    hist(results$results_pooled$record$sigma, xlab="sigma", breaks=15, main="All Clusters")
    hist(results$results_pooled$record$pi, xlab="pi", breaks=15, main="All Clusters")
    hist(results$results_pooled$record$w, xlab="w", breaks=15, main="All Clusters")
    dev.off()}
  }
  
  

 ## Violin plots of estimates
  names_formatted <- substring(names, 9)
  num <- length(results$results[[1]]$record$mu)
  pis <- unlist(lapply(results$results[which.names], function(list) { return(list$record$pi) }))
  
  res_dat <- tibble(cluster = rep(names_formatted, each=num), sample = rep(1:num, length(names)),
                    mu = unlist(lapply(results$results[which.names], function(list) { return(list$record$mu) })),
                    sigma = unlist(lapply(results$results[which.names], function(list) { return(list$record$sigma) })),
                    pi = pis)
  if(coprim.transm==TRUE){
    res_dat$w <- unlist(lapply(results$results[which.names], function(list) { return(list$record$w) }))
  }
  
  # Add pooled trees violin if available
  if (pool.trees==TRUE){
    nump <- nrow(results$results_pooled$record)
    pis <- results$results_pooled$record$pi
    
    temp <- tibble(cluster = rep("ALL", nump), sample = 1:nump,
                   mu = results$results_pooled$record$mu,
                   sigma = results$results_pooled$record$sigma,
                   pi = pis)
    if(coprim.transm==TRUE){
      temp$w <- results$results_pooled$record$w
    }
    res_dat <- rbind(res_dat, temp)
    rm(temp)
  }
  
  getPalette = colorRampPalette(brewer.pal(12, "Set3"))
  res_dat$cluster <- factor(res_dat$cluster, levels = c(mixedsort(names_formatted), "ALL"))
  
  p1 <- ggplot(res_dat, aes(x=cluster, y=mu, fill=cluster)) +
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+
    labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu))
  p1 <- p1 + guides(fill="none") + theme_classic() + 
    scale_fill_manual(values = getPalette(length(names)+1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p2 <- ggplot(res_dat, aes(x=cluster, y=sigma, fill=cluster)) +
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+
    labs(title=paste0(measure, " SD"),x="Cluster", y = expression(sigma))
  p2 <- p2  + guides(fill="none") + theme_classic() + 
    scale_fill_manual(values = getPalette(length(names)+1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p3 <- ggplot(res_dat, aes(x=cluster, y=pi, fill=cluster)) +
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+ylim(0,1)+
    labs(title="Sampling Proportion",x="Cluster", y = expression(pi))
  p3 <- p3  + guides(fill="none") + theme_classic() + 
    scale_fill_manual(values = getPalette(length(names)+1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(coprim.transm==TRUE){
    p6 <- ggplot(res_dat, aes(x=cluster, y=w, fill=cluster)) +
      geom_violin(trim=TRUE)+
      geom_boxplot(width=0.1, fill="white")+ylim(0,1)+
      labs(title="Proportion non-coprimary",x="Cluster", y = "w")
    p6 <- p6  + guides(fill="none") + theme_classic() + 
      scale_fill_manual(values = getPalette(length(names)+1)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  # Print to pdf
  pdf(file = paste0("../Figures/", model, "/violins.pdf"), height = 8.3, width= 11.7, paper="a4r")
  if(coprim.transm==TRUE){grid.arrange(p1, p2, p3, p6, nrow=2, top="Distribution of mean estimates per cluster")
  }else{grid.arrange(p1, p2, p3, nrow=2, top="Distribution of mean estimates per cluster")}
  dev.off()
  if(coprim.transm==TRUE){grid.arrange(p1, p2, p3, p6, nrow=2, top="Distribution of mean estimates per cluster")
  }else{grid.arrange(p1, p2, p3, nrow=2, top="Distribution of mean estimates per cluster")}
  
  
 ## Violin plots of pooled estimates (only)
  if (pool.trees==TRUE){
  names_formatted <- substring(names, 9)
  num <- nrow(results$results_pooled$record)
  pis <- results$results_pooled$record$pi
  res_dat <- tibble(sample = 1:num, cluster = c(rep("ALL", num)),
                     mu = results$results_pooled$record$mu,
                     sigma = results$results_pooled$record$sigma,
                     pi = pis)
  if(coprim.transm==TRUE){
    res_dat$w <- results$results_pooled$record$w
  }
  
  p1 <- ggplot(res_dat, aes(x=cluster, y=mu, fill=cluster)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+
    labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu))
  p1 <- p1 + guides(fill="none") + theme_classic() + scale_fill_brewer(palette="Set3")
  
  p2 <- ggplot(res_dat, aes(x=cluster, y=sigma, fill=cluster)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+
    labs(title=paste0(measure, " SD"),x="Cluster", y = expression(sigma))
  p2 <- p2  + guides(fill="none") + theme_classic() + scale_fill_brewer(palette="Set3")
  
  p3 <- ggplot(res_dat, aes(x=cluster, y=pi, fill=cluster)) + 
    geom_violin(trim=TRUE)+
    geom_boxplot(width=0.1, fill="white")+ylim(0,1)+
    labs(title="Sampling Proportion",x="Cluster", y = expression(pi))
  p3 <- p3  + guides(fill="none") + theme_classic() + scale_fill_brewer(palette="Set3")
  
  if(coprim.transm==TRUE){
    p6 <- ggplot(res_dat, aes(x=cluster, y=w, fill=cluster)) +
      geom_violin(trim=TRUE)+
      geom_boxplot(width=0.1, fill="white")+ylim(0,1)+
      labs(title="Proportion non-coprimary",x="Cluster", y = "w")
    p6 <- p6  + guides(fill="none") + theme_classic() + scale_fill_brewer(palette="Set3")
  }
  
  # Print to pdf
  pdf(file = paste0("../Figures/", model, "/violins_pooled.pdf"), height = 8.3, width= 11.7, paper="a4r")
  if(coprim.transm==TRUE){grid.arrange(p1, p2, p3, p6, nrow=2, top="Distribution of mean estimates per cluster")
  }else{grid.arrange(p1, p2, p3, nrow=2, top="Distribution of mean estimates per cluster")}
  dev.off()
  if(coprim.transm==TRUE){grid.arrange(p1, p2, p3, p6, nrow=2, top="Distribution of mean estimates per cluster")
  }else{grid.arrange(p1, p2, p3, nrow=2, top="Distribution of mean estimates per cluster")}
  }
  
  
  
  
 ## Plot the mean and CI for each cluster - using errorbars

    # Of pooled trees
    if (pool.trees==TRUE){
      res_dat <- tibble(cluster = rep("ALL", 1),
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
        res_dat$w = results$results_pooled$estimates$par[4]
        res_dat$w_ciL = results$results_pooled$estimates$confint[4,1]
        res_dat$w_ciU = results$results_pooled$estimates$confint[4,2]
      }
      
      ncols <- length(names)
      mycols <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
      
      p1 <- ggplot(res_dat, aes(x=cluster, y=mu, col=cluster)) + 
        geom_errorbar(aes(ymin = mu_ciL, ymax = mu_ciU), width=0.2)+
        geom_point()+
        labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu))
      p1 <- p1 + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
      
      p2 <- ggplot(res_dat, aes(x=cluster, y=sigma, col=cluster)) + 
        geom_errorbar(aes(ymin = sigma_ciL, ymax = sigma_ciU), width=0.2)+
        geom_point()
        labs(title=paste0(measure, " SD"),x="Cluster", y = expression(sigma))
      p2 <- p2  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
      
      p3 <- ggplot(res_dat, aes(x=cluster, y=pi, col=cluster)) + 
        geom_errorbar(aes(ymin = pi_ciL, ymax = pi_ciU), width=0.2)+
        geom_point()+ ylim(0,1)+
        labs(title="Sampling Proportion",x="Cluster", y = expression(pi))
      p3 <- p3  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
      
      if (coprim.transm == TRUE){
      p6 <- ggplot(res_dat, aes(x=cluster, y=w, col=cluster)) + 
        geom_errorbar(aes(ymin = w_ciL, ymax = w_ciU), width=0.2)+
        geom_point()+ylim(0,1)+
        labs(title="Proportion non-coprimary",x="Cluster", y = "w")
      p6 <- p6  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
      }
      
      # Print to pdf
      pdf(file = paste0("../Figures/", model, "/bars_pooled.pdf"), height = 8.3, width= 11.7, paper="a4r")
      if (coprim.transm == FALSE){
      grid.arrange(p1, p2, p3, nrow=2, top="Mean and 95% CI parameter estimates (pooled clusters)")
      }else{grid.arrange(p1, p2, p3, p6, nrow=2, top="Mean and 95% CI parameter estimates (pooled clusters)")}
      dev.off()
      if (coprim.transm == FALSE){
        grid.arrange(p1, p2, p3, nrow=2, top="Mean and 95% CI parameter estimates (pooled clusters)")
      }else{grid.arrange(p1, p2, p3, p6, nrow=2, top="Mean and 95% CI parameter estimates (pooled clusters)")}
      
      saveRDS(res_dat, file = paste0("bardata_", model, ".rds"))
      # keep pooled version in memory under new name (for next figure below)
      res_dat_pooled <- res_dat  
      
    }# end of pooled trees section
  
  # Now in each cluster
  res_dat <- tibble(cluster = names_formatted,
                    mu = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$par[1]) })),
                    mu_ciL = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[1,1])})),
                    mu_ciU = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[1,2])})),
                    sigma = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$par[2]) })),
                    sigma_ciL = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[2,1])})),
                    sigma_ciU = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[2,2])})),
                    pi = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$par[3]) })),
                    pi_ciL = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[3,1])})),
                    pi_ciU = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[3,2])})),
                    )
  if (coprim.transm == TRUE){
    res_dat$w = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$par[4]) }))
    res_dat$w_ciL = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[4,1])}))
    res_dat$w_ciU = unlist(lapply(results$results[which.names], function(list) { return(list$estimates$confint[4,2])}))
  }
  
  # Add pooled trees if available
  ncols <- length(names) 
  if (pool.trees==TRUE){
  ncols <- length(names) + 1
  res_dat <- rbind(res_dat, res_dat_pooled)
  }
  mycols <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
  res_dat$cluster <- factor(res_dat$cluster, levels = c(mixedsort(names_formatted), "ALL"))
  
  p1 <- ggplot(res_dat, aes(x=cluster, y=mu, col=cluster)) +
    geom_errorbar(aes(ymin = mu_ciL, ymax = mu_ciU), width=0.2)+
    geom_point()+ ylim(0,max(11, max(res_dat$mu_ciU))) +
    labs(title=paste0(measure, " Mean"),x="Cluster", y = expression(mu)) 
  p1 <- p1 + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14))
  
  
  p2 <- ggplot(res_dat, aes(x=cluster, y=sigma, col=cluster)) +
    geom_errorbar(aes(ymin = sigma_ciL, ymax = sigma_ciU), width=0.2)+
    geom_point()+ ylim(0,max(7, max(res_dat$sigma_ciU))) +
    labs(title=paste0(measure, " Standard Deviation"),x="Cluster", y = expression(sigma))
  p2 <- p2  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14))
  
  p3 <- ggplot(res_dat, aes(x=cluster, y=pi, col=cluster)) +
    geom_errorbar(aes(ymin = pi_ciL, ymax = pi_ciU), width=0.2)+
    geom_point()+ylim(0,1)+
    labs(title="Sampling Proportion",x="Cluster", y = expression(pi))
  p3 <- p3  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14))
  
  if (coprim.transm == TRUE){
    p6 <- ggplot(res_dat, aes(x=cluster, y=w, col=cluster)) + 
      geom_errorbar(aes(ymin = w_ciL, ymax = w_ciU), width=0.2)+
      geom_point()+ylim(0,1)+
      labs(title="Proportion Non-coprimary",x="Cluster", y = "w")
    p6 <- p6  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14))
  }
  
  # Print to pdf
  pdf(file = paste0("../Figures/", model, "/bars.pdf"),  height = 8.3, width= 11.7, paper="a4r")
  if (coprim.transm == FALSE){
    grid.arrange(p1, p2, p3, nrow=2, top= textGrob(paste0("Wave ", which.wave), x = 0.13, gp=gpar(fontsize=14, fontface = "bold")))
  }else{grid.arrange(p1, p2, p3, p6, nrow=2, top= textGrob(paste0("Wave ", which.wave), x = 0.07, gp=gpar(fontsize=16, fontface = "bold")))}
  dev.off()
  if (coprim.transm == FALSE){
    grid.arrange(p1, p2, p3, nrow=2, top= textGrob(paste0("Wave ", which.wave), x = 0.13, gp=gpar(fontsize=14, fontface = "bold")))
  }else{grid.arrange(p1, p2, p3, p6, nrow=2, top= textGrob(paste0("Wave ", which.wave), x = 0.07, gp=gpar(fontsize=16, fontface = "bold")))}
  
   
  ## A similar plot, but for a selection of trees from a single cluster
  z <- qnorm(0.95+(1-0.95)/2)
  # Sample 10 of the trees at random
  samp <- sample(1:length(results$results[[1]]$record$mu), 10)
  res_dat2 <- tibble(sample = as.factor(c(1:10)),
                     mu = results$results[[which.names[1]]]$record$mu[samp],
                     mu_ciL = results$results[[which.names[1]]]$record$mu[samp]-z*results$results[[which.names[1]]]$record$se.mu[samp],
                     mu_ciU = results$results[[which.names[1]]]$record$mu[samp]+z*results$results[[which.names[1]]]$record$se.mu[samp],
                     sigma = results$results[[which.names[1]]]$record$sigma[samp],
                     sigma_ciL = results$results[[which.names[1]]]$record$sigma[samp]-z*results$results[[which.names[1]]]$record$se.sigma[samp],
                     sigma_ciU = results$results[[which.names[1]]]$record$sigma[samp]+z*results$results[[which.names[1]]]$record$se.sigma[samp],
                     pi = results$results[[which.names[1]]]$record$pi[samp],
                     pi_ciL = results$results[[which.names[1]]]$record$pi[samp]-z*results$results[[which.names[1]]]$record$se.pi[samp],
                     pi_ciU = results$results[[which.names[1]]]$record$pi[samp]+z*results$results[[which.names[1]]]$record$se.pi[samp]
  )
  if (coprim.transm == TRUE){
  res_dat2$w = results$results[[which.names[1]]]$record$w[samp]
  res_dat2$w_ciL = results$results[[which.names[1]]]$record$w[samp]-z*results$results[[which.names[1]]]$record$se.w[samp]
  res_dat2$w_ciU = results$results[[which.names[1]]]$record$w[samp]+z*results$results[[which.names[1]]]$record$se.w[samp]
  }
  
  mycols <- colorRampPalette(brewer.pal(8, "Set2"))(10)
  
  p1 <- ggplot(res_dat2, aes(x=sample, y=mu, col=sample)) +
    geom_errorbar(aes(ymin = mu_ciL, ymax = mu_ciU), width=0.2)+
    geom_point()+
    labs(title=paste0(measure, " Mean"),x="Sampled tree", y = expression(mu))
  p1 <- p1 + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
  
  p2 <- ggplot(res_dat2, aes(x=sample, y=sigma, col=sample)) +
    geom_errorbar(aes(ymin = sigma_ciL, ymax = sigma_ciU), width=0.2)+
    geom_point()+
    labs(title=paste0(measure, " SD"),x="Sampled tree", y = expression(sigma))
  p2 <- p2  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
  
  p3 <- ggplot(res_dat2, aes(x=sample, y=pi, col=sample)) +
    geom_errorbar(aes(ymin = pi_ciL, ymax = pi_ciU), width=0.2)+
    geom_point()+ylim(0,1)+
    labs(title="Sampling Proportion",x="Sampled tree", y = expression(pi))
  p3 <- p3  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
  
  p6 <- ggplot(res_dat2, aes(x=sample, y=w, col=sample)) + 
    geom_errorbar(aes(ymin = w_ciL, ymax = w_ciU), width=0.2)+
    geom_point()+ylim(0,1)+
    labs(title="Proportion non-coprimary",x="Sampled tree", y = "w")
  p6 <- p6  + guides(col="none") + theme_classic() + scale_colour_manual(values = mycols)
  
  
  # Print to pdf
  pdf(file = paste0("../Figures/", model, "/bars_comparetrees.pdf"), height = 8.3, width= 11.7, paper="a4r")
  if (coprim.transm == FALSE){
    grid.arrange(p1, p2, p3, nrow=2, top="Mean and 95% CI parameter estimates from 10 random sampled trees")
  }else{grid.arrange(p1, p2, p3, p6, nrow=2, top="Mean and 95% CI parameter estimates from 10 random sampled trees")}
  dev.off()
  if (coprim.transm == FALSE){
    grid.arrange(p1, p2, p3, nrow=2, top="Mean and 95% CI parameter estimates from 10 random sampled trees")
  }else{grid.arrange(p1, p2, p3, p6, nrow=2, top="Mean and 95% CI parameter estimates from 10 random sampled trees")}
  
  
  
  ## If using a coprimary model, what is the overall sampling probabilty?
  if (coprim.transm==TRUE){
    ## Sampling frac in the coprimary model 
    # w*pi + (1-w)*2/3
    p_overall <- 0
    for (i in 1:length(names)){
      p_overall[i] <- results$results[[which.names[i]]]$estimates$par[4]*results$results[[which.names[i]]]$estimates$par[3] + (1-results$results[[which.names[i]]]$estimates$par[4])*2/3
    }
    print("Overall sampling probability, w*pi + (1-w)*2/3, per cluster is ")
    print(p_overall)
  }
  
  
  
  
  ## Plot densities of the estimates (of a single SI per cluster)
  ncols <- length(names) + 1
  mycols <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
  mycols[length(names) + 1] <- c("Black")
  
  tmax <- 0
  for (i in 1:length(names)){
    tmax <- max(tmax,results$results[[which.names[i]]]$estimates$par[1]+ 6*results$results[[which.names[i]]]$estimates$se[1])
  }
  dat <- tibble(time = seq(0,tmax, length.out=200))
  for (i in 1:length(names)){
    dat[[names[i]]] <- dgamma(dat$time, 
                              shape = results$results[[which.names[i]]]$estimates$par[1]^2/results$results[[which.names[i]]]$estimates$par[2]^2,
                              rate = results$results[[which.names[i]]]$estimates$par[1]/results$results[[which.names[i]]]$estimates$par[2]^2)
  }
  
  # Add pooled analysis where available
  if (pool.trees==TRUE){
    dat[["ALL"]] <- dgamma(dat$time,
                           shape = results$results_pooled$estimates$par[1]^2/results$results_pooled$estimates$par[2]^2,
                           rate = results$results_pooled$estimates$par[1]/results$results_pooled$estimates$par[2]^2)
  }
  
  # change to long format              
  datl <- dat %>% pivot_longer(cols = c(2:ncol(dat)))
  datl <- datl %>% rename(Cluster = name, density = value)
  
  if (pool.trees==TRUE){
    datl$Cluster <- factor(datl$Cluster, levels=c(mixedsort(names),
                                                  "ALL"))}
  
  #Larger line size for pooled trees
  l <- ifelse(datl$Cluster=="ALL", 2, 1)
  datl$size <- l
  # Legend labels
  if(pool.trees==TRUE){leglabs <- c(str_to_title(gsub("_", " ", mixedsort(names), fixed=TRUE)), "ALL")
  }else{leglabs <-str_to_title(gsub("_", " ", mixedsort(names), fixed=TRUE))}
  
  ps <- ggplot(datl, aes(time, density, group=Cluster, color = Cluster, size = size)) + 
    geom_line() + 
    theme_minimal() +
    ylab("Density") + xlab("Time (days)") +
    scale_color_manual(values = mycols, labels = leglabs) + 
    scale_size(range = c(0.5, 2), guide="none") + ggtitle(paste0("Wave ", which.wave)) +
    theme(text = element_text(size=24), plot.title = element_text(size = 26, face = "bold"))
  
  ggsave(paste0("../Figures/", model, "/SI_densities.pdf"),ps, width=11.7, height=8.3, units="in", paper="a4r")
  
  ps  
  
  return(list(Table = Table))
}
  