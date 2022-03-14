#############################################
### Plot epi curves of Melbourne clusters ###
#############################################

epi_curves <- function(metadata, which.wave=1){
  
  
  # Desired cluster names
  names <- as.character(mixedsort(unique(metadata$cluster_id)))

  # Get number of cases per day by symptom onset, and add cumulative sum column
  g_data <- metadata %>% filter(cluster_id %in% names) %>%
            group_by(cluster_id, onset_date) %>% summarise(n = n()) %>% 
            #count(onset_date, cluster_id)
            group_by(cluster_id) %>%
            mutate(csum = cumsum(n))
  g_data

  # Plot
  g_data$cluster_id <- factor(as.character(g_data$cluster_id), levels=mixedsort(names))
  g_data$onset_date  <- as.Date(g_data$onset_date )
  
  getPalette <- colorRampPalette(brewer.pal(8, "Set1"))(length(names)+1)
  
  p1 <- g_data %>%
   ggplot( aes(x=onset_date, y=csum, group=cluster_id, color=cluster_id)) +
   geom_line() +
    theme_minimal() +
    ggtitle(paste0("Wave ", which.wave)) +
    ylab("Number of cases (cumulative)") + xlab("Symptom onset date") + 
    guides(color=guide_legend(title="Cluster")) +
    scale_colour_manual(values = getPalette) + theme(text = element_text(size=12))
  
  # Print to pdf
  pdf(paste0("../Figures/symptom_epicurves_w", which.wave, ".pdf", height = 6, width= 5))
  p1
  dev.off()


  
return(p1)  
}


