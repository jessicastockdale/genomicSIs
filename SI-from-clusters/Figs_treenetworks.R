#############################################
### Plot tree networks of Melbourne clusters ###
#############################################

treenets <- function(metadata, names, coprim.transm, pi.model, w.model, threshold=0.1, which.wave=1){
  
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
  print(paste0("Confirming, chosen model is: ", model)) 
  # Create directory for the figures, if it doesn't already exist
  dir.create(file.path(paste0("../Figures/Sampledtreenets")), showWarnings = FALSE)
  
  
  ## Read in sampled tree files
  sampled_trees <- vector(mode = "list", length = length(names))
  
    for (i in 1:length(names)){
      sampled_trees[[i]] <- readRDS(file = paste0("sampled_trees/", names[i], model, ".rds"))
    }


## Build a matrix per cluster of rows = infector, columns = infectee, 
#  number  = prop of trees with that pair
M <- vector(mode = "list", length = length(names))
for (i in 1:length(names)){
  # Set up the matrix
  nams <- metadata$case_id[metadata$cluster_id==sub(".*_", "", names[i])]
  M[[i]] <- matrix(0, nrow = length(nams), ncol = length(nams))
  rownames(M[[i]]) <- nams
  colnames(M[[i]]) <- nams
  
  # Fill in with proportion of trees that include each infector/ee pair
  for (j in 1:length(sampled_trees[[i]])){
    for (k in 1:nrow(sampled_trees[[i]][[1]])){
         M[[i]][sampled_trees[[i]][[j]][k,c(2,3)][[1]], sampled_trees[[i]][[j]][k,c(2,3)][[2]]] <-
           M[[i]][sampled_trees[[i]][[j]][k,c(2,3)][[1]], sampled_trees[[i]][[j]][k,c(2,3)][[2]]] +
           1
    }
  }
  M[[i]] <- M[[i]]/length(sampled_trees[[i]])
  # Sums of columns should all be 1, or 0 (those that are always initial infectives)
  # Each i should return TRUE
  print(all(colSums(M[[i]])==1 | colSums(M[[i]])==0))
}

## Set all entries below the threshold to 0
M2 <- vector(mode = "list", length = length(names))
for (i in 1:length(names)){
  M2[[i]] <- M[[i]]
  M2[[i]][M2[[i]] < threshold] <- 0
}

## Plot as networks
pl <- vector(mode = "list", length = length(names))
plot_title <- str_to_title(gsub("_", " ", names, fixed=TRUE))
for (i in 1:length(names)){
  ## make igraph object
  my.gr <- graph_from_adjacency_matrix(M2[[i]], mode="directed", weighted=T)
  
  ## convert to VisNetwork-list
  my.visn <- toVisNetworkData(my.gr)
  ## copy column "weight" to new column "value" in list "edges"
  my.visn$edges$value <- my.visn$edges$weight
  
  my.visn$edges
  
  pl[[i]] <- visNetwork(my.visn$nodes, my.visn$edges, main=plot_title[i]) %>%
    visEdges(arrows =list(to = list(enabled = TRUE, scaleFactor = 1))) %>%
    visLayout(randomSeed = 12) %>%
    visIgraphLayout(layout = "layout_nicely", type = "full")

  visSave(pl[[i]], file = paste0("../Figures/Sampledtreenets/", model, "_",names[i], ".html"))
}


return(pl)  
}



