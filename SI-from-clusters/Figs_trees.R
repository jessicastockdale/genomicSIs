#######################################
### Plots of some our sampled trees ###
#######################################

samptrees <- function(metadata, names, coprim.transm, pi.model, w.model, which.wave=1){
  
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
  if (which.wave==2){
    model <- paste0(model, "_w2")
  }
  measure <- "Serial Interval"
  print(paste0("Confirming, chosen model is: ", model)) 
  
  ## Read in sampled tree files
  sampled_trees <- vector(mode = "list", length = length(names))
    for (i in 1:length(names)){
      sampled_trees[[i]] <- readRDS(file = paste0("sampled_trees/", names[i], model, ".rds"))
    }
  
  
  # Create output folder, if it doesn't already exist
  dir.create("../Figures/Samp_trees", showWarnings = FALSE)
  # Create colour palette
  ncols <- length(names)
  mycols <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)

for (i in 1:length(sampled_trees)){
  for (j in 1:length(sampled_trees[[i]])){
    which_tree <- sampled_trees[[i]][[j]]  # loops over all trees

    edges <- data.frame(from = which_tree$case_i,to =which_tree$case_j ,weight = which_tree$onset_diff)
    nodes <- data.frame(label= unique(c(which_tree$case_i, which_tree$case_j)))

    routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
    if (which.wave==1){
      pdf(file = paste0("../Figures/Samp_trees/", names[i], "_t", j, ".pdf"), width = 8, height = 8) 
    }else{pdf(file = paste0("../Figures/Samp_trees/w2_", names[i], "_t", j, ".pdf"), width = 8, height = 8) }
      par(mar=c(0,0,0,0)+.1)
    plot(routes_igraph, layout = layout_with_graphopt, edge.arrow.size = 0.5, 
     vertex.label.cex=1.2, vertex.size=10, vertex.color = mycols[i],
     vertex.frame.color = mycols[i],
     edge.color=mycols[i], vertex.label.color="black", edge.label = edges$weight,
     edge.label.color="darkgray", edge.label.cex = 1.5, edge.width = 1.2)
    title(paste0(str_to_title(gsub("_", " ", names[i], fixed=TRUE)), " Tree ", j), adj=0, line=-1, cex.main=1.5)
    dev.off()

    }
  }

}

