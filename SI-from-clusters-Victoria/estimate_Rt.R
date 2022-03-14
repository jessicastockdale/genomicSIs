# Cluster-specific Rt estimates

estimate.rt <- function(metadata, coprim.transm,  pi.model="none", w.model="none", pool.trees = FALSE, which.wave=1){
  
  
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
print(paste0("Confirming, chosen model is: ", model)) 


## Results: SERIAL INTERVALS BY CLUSTER
SI_res <- readRDS(file = paste0("resultstable_", model, ".rds"))
# Don't run 'All' clusters result, if it's present
if ("ALL" %in% names(SI_res)){
  SI_res <- SI_res[,-ncol(SI_res)]
}

## Metadata: INCIDENCE BY CLUSTER
cases <- metadata %>%
 filter(cluster_id %in% sub(".*?_", "", names(SI_res))) %>%
 filter(!is.na(onset_date)) %>%
 group_by(onset_date, cluster_id) %>%
 summarize(incidence=n()) 

cases$onset_date <- as.Date(cases$onset_date)
cases$cluster_id <- factor(cases$cluster_id)
cases %>%
  ggplot(aes(x=onset_date, y=incidence, group=cluster_id, fill=cluster_id))+
  geom_col()+
  theme_classic()


ncols <- ncol(data.frame(SI_res))+1
mycols <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)


## Run: OVER ALL CLUSTERS

get_Rt <- function(cl, with.inset=FALSE, compare=FALSE){
  # compare=TRUE will plot a second Rt estimate, using published serial interval parameters
  # Using Bi et al 2020, gamma mean 6?3 days, SD 4?2 days. https://doi.org/10.1016/S1473-3099(20)30287-5
  
  
  # Note: can also take si_distr, min/max
  config <- make_config(list(mean_si= SI_res['mu', cl], 
                            std_si = SI_res['sigma', cl]))
  
  if (compare==TRUE){
    config2 <- make_config(list(mean_si= 6.3, 
                               std_si = 4.2))
  }

  # Run Rt estimation
  data <- filter(cases, cluster_id == sub(".*?_", "", cl))
  
  #Create incidence data frame in required format
  incid <- data.frame(I = data$incidence, dates = data$onset_date)
  incid <- incid %>% pad()
  incid$I[is.na(incid$I)] <- 0
  
  resR <- estimate_R(incid,
    method='parametric_si',
    config=config)
  
  if (compare==TRUE){
    resR2 <- estimate_R(incid,
                       method='parametric_si',
                       config=config2)
  }
  
  #plot
  if (compare==FALSE){
  gg <- plot(resR, "R", legend=FALSE, options_R = list(col = mycols[cl==names(SI_res)])) + 
    ggtitle(paste0('Cluster ', sub(".*?_", "", cl), '. Size ', sum(data$incidence))) +
    labs(y="R", x="")+theme_classic()+
    theme(legend.position='none',
    text = element_text(size=14),
    axis.text.y=element_text(size=14),
    axis.text.x=element_text(size=14, angle = 45, hjust = 1))
  if (with.inset){
    inset <- plot(resR, "incid", legend=FALSE, options_I = list(col = c("grey49"))) +theme_classic()+ scale_y_continuous(breaks= pretty_breaks())+
      theme(legend.position='none',
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=10), plot.margin=grid::unit(c(0,0,0,0), "mm")) +
      ggtitle("")+labs(x="", y="")
    gg <- ggdraw() + draw_plot(gg) + draw_plot(inset, x=0.65, y=0.6, width=0.3, height=0.3)
  }
  }else{
    gg <- estimate_R_plots(list(resR2, resR), "R", legend=TRUE, options_R = list(col = c("darkgray",mycols[cl==names(SI_res)][1]))) + 
      ggtitle(paste0('Cluster ', sub(".*?_", "", cl), '. Size ', sum(data$incidence))) +
      labs(y="R", x="")+theme_classic()+
      theme(legend.position='none',
            text = element_text(size=14),
            axis.text.y=element_text(size=14),
            axis.text.x=element_text(size=14, angle = 45, hjust = 1)) 
    if (with.inset){
      inset <- plot(resR, "incid", legend=FALSE, options_I = list(col = c("grey49"))) +theme_classic()+
        theme(legend.position='none',
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.y=element_text(size=10), plot.margin=grid::unit(c(0,0,0,0), "mm")) +
        ggtitle("")+labs(x="", y="")
      gg <- ggdraw() + draw_plot(gg) + draw_plot(inset, x=0.65, y=0.6, width=0.3, height=0.3)
    }
    
  }
  return(gg)
}

get_incid <- function(cl, with.inset=FALSE){
 
  # Note: can also take si_distr, min/max
  config <- make_config(list(mean_si= SI_res['mu', cl], 
                            std_si = SI_res['sigma', cl]))
  # Run Rt estimation
  data <- filter(cases, cluster_id == sub(".*?_", "", cl))

  #Create incidence data frame in required format
  incid <- data.frame(I = data$incidence, dates = data$onset_date)
  incid <- incid %>% pad()
  incid$I[is.na(incid$I)] <- 0
  
  resR <- estimate_R(incid,
                    method='parametric_si',
                    config=config)

  #plot
  gg <- plot(resR, "incid", legend=FALSE, options_I = list(col = mycols[cl==names(SI_res)])) + 
    ggtitle(paste0('Cluster ', sub(".*?_", "", cl), '. Size ', sum(data$incidence))) +
    labs(y="Cases", x="")+theme_classic()+
    scale_x_date(breaks = "week", labels=date_format("%m-%d"))+
    theme(legend.position='none',
    text = element_text(size=10),
    axis.text.y=element_text(size=8),
    axis.text.x=element_text(size=8, angle = 45, hjust = 1)) 
  return(gg)
}

get_both <- function(cl, compare=FALSE){
    # compare=TRUE will plot a second Rt estimate, using published serial interval parameters
    # Using Bi et al 2020, gamma mean 6.3 days, SD 4.2 days. https://doi.org/10.1016/S1473-3099(20)30287-5
    
    config <- make_config(list(mean_si= SI_res['mu', cl], 
                               std_si = SI_res['sigma', cl]))
    if (compare==TRUE){
      config2 <- make_config(list(mean_si= 6.3, 
                                  std_si = 4.2))
    }
    
    # Run Rt estimation
    data <- filter(cases, cluster_id == sub(".*?_", "", cl))
    
    #Create incidence data frame in required format
    incid <- data.frame(I = data$incidence, dates = data$onset_date)
    incid <- incid %>% pad()
    incid$I[is.na(incid$I)] <- 0
    
    resR <- estimate_R(incid,
                       method='parametric_si',
                       config=config)
    
    if (compare==TRUE){
      resR2 <- estimate_R(incid,
                          method='parametric_si',
                          config=config2)
    }

    
    #plot
    df_incid <- tibble(Date = resR$dates, Incidence = resR$I_local + resR$I_imported)
    if (compare==FALSE){
      gg <- plot(resR, "R", legend=FALSE, options_R = list(col = mycols[cl==names(SI_res)])) + 
        ggtitle(paste0('Cluster ', sub(".*?_", "", cl), '. Size ', sum(data$incidence))) +
        labs(y="R", x="")+theme_classic()+
        theme(legend.position='none',
              text = element_text(size=14),
              axis.text.y=element_text(size=14),
              axis.text.x=element_text(size=14, angle = 45, hjust = 1)) +
        geom_col(data=df_incid, aes(x=Date, y=Incidence), alpha=0.1) + scale_y_continuous(breaks= pretty_breaks(), sec.axis = dup_axis(name = "Incidence"))
    }else{
      gg <- estimate_R_plots(list(resR2, resR), "R", legend=TRUE, options_R = list(col = c("skyblue2",mycols[cl==names(SI_res)][1]))) + 
        ggtitle(paste0('Cluster ', sub(".*?_", "", cl), '. Size ', sum(data$incidence))) +
        labs(y="Rt", x="")+theme_classic()+
        theme(legend.position='none',
              text = element_text(size=14),
              axis.text.y=element_text(size=14),
              axis.text.x=element_text(size=14, angle = 45, hjust = 1))+
        geom_col(data=df_incid, aes(x=Date, y=Incidence), alpha=0.1) + scale_y_continuous(breaks= pretty_breaks(), sec.axis = dup_axis(name = "Incidence"))
      
    }
    
  
  return(gg)
}

get_Rt_point <- function(cl, t=5, compare=FALSE){ 
  # Note: can also take si_distr, min/max
  config <- make_config(list(mean_si= SI_res['mu', cl], 
                            std_si = SI_res['sigma', cl]))
  
  
  if (compare==TRUE){
    config2 <- make_config(list(mean_si= 6.3, 
                                std_si = 4.2))
  }
  
  # Run Rt estimation
 # if (cl!="All"){
    data <- filter(cases, cluster_id == sub(".*?_", "", cl))
 # }else{ # collapse down over clusters
 #   data <- cases
 #   data <- data %>%
 #     group_by(onset_date) %>%
 #     summarize(incidence=sum(incidence)) }
  
  #Create incidence data frame in required format
  incid <- data.frame(I = data$incidence, dates = data$onset_date)
  incid <- incid %>% pad()
  incid$I[is.na(incid$I)] <- 0
  
  resR <- estimate_R(incid,
                    method='parametric_si',
                    config=config)
  if (compare==TRUE){
    resR2 <- estimate_R(incid,
                        method='parametric_si',
                        config=config2)
  }
  
  if (t=="max.inc"){
    maxinc <- which.max(incid$I)
    R <- resR$R[maxinc,]
    R$Method <- c("Fixed")
    if (compare==TRUE){R[2,] <- resR2$R[maxinc,]
                       R$Method <- c("Estimated", "Fixed")}
  }else{R <- resR$R[t,]
  R$Method <- c("Fixed")
  if (compare==TRUE){R[2,] <- resR2$R[t,]
                     R$Method <- c("Estimated", "Fixed")}
  }
  
  R$Cluster <- sub(".*?_", "", cl)
  return(R)
}


# Generate various different versions of the Rt figure to choose between:

resR <- lapply(names(SI_res), function(p) {
  tryCatch(get_Rt(p, compare=TRUE, with.inset=FALSE), error = function(e) return(NULL))})
gg_res <- ggpubr::ggarrange(plotlist=resR, ncol=4, nrow=3)
pdf(paste0("../Figures/", model, "/Rt_ests.pdf"), width=11.7, height=8.3, paper="a4r")
print(gg_res)
dev.off()

resRi <- lapply(names(SI_res), function(p) {
  tryCatch(get_Rt(p, compare=TRUE, with.inset=TRUE), error = function(e) return(NULL))})
gg_resi <- ggpubr::ggarrange(plotlist=resRi, ncol=4, nrow=3)
pdf(paste0("../Figures/", model, "/Rt_ests_inset.pdf"), width=11.7, height=8.3, paper="a4r")
print(gg_resi)
dev.off()

resR2 <- lapply(names(SI_res), function(p) {
  tryCatch(get_incid(p), error = function(e) return(NULL))} )
gg_res2 <- ggpubr::ggarrange(plotlist=resR2, ncol=4, nrow=3)
pdf(paste0("../Figures/", model, "/cluster_incidence.pdf"), width=11.7, height=8.3, paper="a4r")
print(gg_res2)
dev.off()

resR3 <- lapply(names(SI_res), function(p) {
  tryCatch(get_both(p, compare=TRUE), error = function(e) return(NULL))} )
gg_res3 <- ggpubr::ggarrange(plotlist=resR3, ncol=4, nrow=3)
pdf(paste0("../Figures/", model, "/Rt_ests_with_incidence.pdf"), width=11.7, height=8.3, paper="a4r")
print(gg_res3)
dev.off()

## Rt's at one time point - between 6 and 12 days
resR4 <- names(SI_res) %>% map_dfr(function(p) {
  tryCatch(get_Rt_point(p , compare=TRUE), error = function(e) return(NULL))})
names(resR4) <- str_remove(names(resR4), "\\(R\\)")
# Plot point estimates (moving avg. between days 6 and 12)
gg_res4 <- ggplot(resR4, aes_string(x='Cluster', y='Median', col='Cluster', group="Method"))+
  geom_point(aes_string(shape="Method"), size=2.5, position=position_dodge(width=0.8))+
  scale_shape_manual(values=c(16,1)) +
  geom_errorbar(aes_string(ymin="Quantile.0.05", ymax="Quantile.0.95", linetype="Method"), size=0.6, 
                position=position_dodge(width=0.8))+
   labs(y="Rt", x="Cluster", title="Rt Estimate between 6 and 12 days")+
    theme_classic()+ scale_colour_manual(values = mycols)+
    theme(text = element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.text.x=element_text(size=12, angle=45, vjust = 0.25, hjust=0.5))+
    guides(col=guide_legend(ncol=1))
gg_res4
ggsave(paste0("../Figures/", model, "/compare_Rts.pdf"), gg_res4, width=11.7, height=8.3, units="in", paper="a4r")
  

return(list(rt_figs = gg_resi, incidence_fig = gg_res2, rt_comparison = gg_res4))
}
