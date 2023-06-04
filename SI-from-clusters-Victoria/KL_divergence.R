
### Calculate KL divergence of prior/posterior and likelihood/posterior, for our sampled trees

# Giri Gopalan, https://arxiv.org/pdf/1511.01214.pdf
# prior information = DKL(posterior, likelihood)
# likelihood information = DKL(posterior, prior) 

library(LaplacesDemon)

### Set up
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
# Configuration of the SI estimation
init_ <- list(mu=8,sigma=2,pi=0.8) #initial values for alpha, beta, pi
low_ <- c(0.001, 0.001, 0.3) #lower bound of each estimate
upp_ <- c(30, 30, 1.001) #upper bound of each estimate
cores_ <- availableCores() - 1 #how many core processors
# initial, upper and lower values for w where needed
if (coprim.transm==TRUE){init_$w <- 0.8
low_ <- c(low_, 0.01)
upp_ <- c(upp_, 1.0)} 




### KL calculations. Matrices in terms of pi, w


# The beta prior: Beta(12,11)
the_vec <- seq(0,1,0.001)
prior_Y <- dbeta(the_vec, 12, 11)


# Load results (MAPs in each tree, in each cluster)
results_sam1 <- readRDS(file = paste0(model, "_w", which.wave, "_results.rds"))

# Structure for results
lh_info <- matrix(NA, length(names), 10)
prior_info <- matrix(NA, length(names), 10)

# Required functions
findpost_pi <- function(k, pars){
  pars2 <- c(pars[1], pars[2], k, pars[4])
  # Post
  post_X <- postOPTIM(pars2, si, pi_prior=pi.info, w_prior=w.info)
}
findlike_pi <- function(k, pars){
  pars2 <- c(pars[1], pars[2], k, pars[4])
  # Lh
  like_Y <- llOPTIM(pars2, si)
}

findpost_w <- function(k, pars){
  pars2 <- c(pars[1], pars[2], pars[3], k)
  # Post
  post_X <- postOPTIM(pars2, si, pi_prior=pi.info, w_prior=w.info)
}
findlike_w <- function(k, pars){
  pars2 <- c(pars[1], pars[2], pars[3], k)
  # Lh
  like_Y <- llOPTIM(pars2, si)
}


### Run it: PI PRIOR
for (i in 1:length(names)){
 
    print(paste0("Currently running ", names[i]))

    # Load sampled trees
    trees_sam <- readRDS(file = paste0("sampled_trees/",names[i], model, "_w", which.wave, ".rds"))
    
    # Which trees got used?
    which.samples<-rownames(results_sam1[[i]]$record)
    which.samples <- as.numeric(gsub("mu", "", which.samples))
    if(is.na(which.samples[1])){which.samples[1] <- 0}
    
    count <- 1
    
    # Sample 10 of the trees to compute
    for (j in sample(0:(length(trees_sam)-1))){
      
    if (count<11){
      
    print(paste0("Currently running tree", j))
      
    # Did this tree get used?
    if (j %in% which.samples){
    
      # Calculate posterior and likelihood densities
      treedata=trees_sam[[j+1]]
      si <- treedata[,"onset_diff"]
      si <- si[!is.na(si)]
      pars = unlist(init_) # Will change pi/w to calculate the KL. Set remaining params to their maps
      pars[1] = results_sam1[[i]]$record$mu[which(which.samples==j)]
      pars[2] = results_sam1[[i]]$record$sigma[which(which.samples==j)]
      # pars[3] = results_sam1[[i]]$record$pi[which(which.samples==j)] this one is varied here
      pars[4] = results_sam1[[i]]$record$w[which(which.samples==j)]

      post_X <- sapply(the_vec, findpost_pi, pars = pars)
      like_Y <- sapply(the_vec, findlike_pi, pars = pars)
      
      # Normalize
      post_X <- 1000*exp(post_X)/sum(exp(post_X))
      like_Y <- 1000*exp(like_Y)/sum(exp(like_Y))
      
      # Calculate KL divergence
      postprior <- KLD(post_X, prior_Y)
      postlike <- KLD(post_X, like_Y)
      lh_info[i,count] <- postprior$sum.KLD.px.py
      prior_info[i,count] <- postlike$sum.KLD.px.py
      
      count <- count + 1
      
    }
      
    }
      
      plot(the_vec, prior_Y)
      plot(the_vec, like_Y)
      plot(the_vec, post_X)
    
    }

  
}

lh_info_pi <- lh_info
prior_info_pi <- prior_info

### Run it: W PRIOR
for (i in 1:length(names)){
  
  print(paste0("Currently running ", names[i]))
  
  # Load sampled trees
  trees_sam <- readRDS(file = paste0("sampled_trees/",names[i], model, "_w", which.wave, ".rds"))
  
  # Which trees got used?
  which.samples<-rownames(results_sam1[[i]]$record)
  which.samples <- as.numeric(gsub("mu", "", which.samples))
  if(is.na(which.samples[1])){which.samples[1] <- 0}
  
  count <- 1
  
  # Sample 10 of the trees to compute
  for (j in sample(0:(length(trees_sam)-1))){
    
    if (count<11){
      
      print(paste0("Currently running tree", j))
      
      # Did this tree get used?
      if (j %in% which.samples){
        
        # Calculate posterior and likelihood densities
        treedata=trees_sam[[j+1]]
        si <- treedata[,"onset_diff"]
        si <- si[!is.na(si)]
        pars = unlist(init_) # Will change pi/w to calculate the KL. Set remaining params to their maps
        pars[1] = results_sam1[[i]]$record$mu[which(which.samples==j)]
        pars[2] = results_sam1[[i]]$record$sigma[which(which.samples==j)]
        pars[3] = results_sam1[[i]]$record$pi[which(which.samples==j)]  
        #pars[4] = results_sam1[[i]]$record$w[which(which.samples==j)] this one is varied here
        
        post_X <- sapply(the_vec, findpost_w, pars = pars)
        like_Y <- sapply(the_vec, findlike_w, pars = pars)
        
        # Normalize
        post_X <- 1000*exp(post_X)/sum(exp(post_X))
        like_Y <- 1000*exp(like_Y)/sum(exp(like_Y))
        
        # Calculate KL divergence
        postprior <- KLD(post_X, prior_Y)
        postlike <- KLD(post_X, like_Y)
        lh_info[i,count] <- postprior$sum.KLD.px.py
        prior_info[i,count] <- postlike$sum.KLD.px.py
        
        count <- count + 1
        
      }
      
    }
    
    plot(the_vec, prior_Y)
    plot(the_vec, like_Y)
    plot(the_vec, post_X)
    
  }
  
  
}


lh_info_w <- lh_info
prior_info_w <- prior_info




# Plot the results (% stacked bars)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gtools)

new.names <- c("A7",  "A8",  "A1",  "A10", "A5",  "A2",  "A4",  "A3",
               "A9", "A6")



tlh_info_pi <- as_tibble(lh_info_pi)
colnames(tlh_info_pi) <- paste0("Tree ", 1:10)
pilh_res<-pivot_longer(cbind(new.names,tlh_info_pi), cols = 2:11)
pilh_res$Source <- c(rep("Likelihood", nrow(pilh_res)))

tp_info_pi <- as_tibble(prior_info_pi)
colnames(tp_info_pi) <- paste0("Tree ", 1:10)
pip_res<-pivot_longer(cbind(new.names,tp_info_pi), cols = 2:11)
pip_res$Source <- c(rep("Prior", nrow(pip_res)))

pi_res <- rbind(pilh_res, pip_res)

pi_res$new.names <- factor(pi_res$new.names, levels = mixedsort(new.names))

# Stacked + percent. Average over the 10 trees
pi_res <- pi_res %>% group_by(new.names, name) %>% pivot_wider(names_from = "Source") %>% mutate(Likelihood = Likelihood/(Likelihood+Prior),
                                                                                               Prior = Prior/(Likelihood+Prior)) %>% 
  pivot_longer(cols = c("Likelihood", "Prior"), names_to = "Source")


 p1 <- pi_res %>% group_by(new.names, Source) %>% summarize(value = mean(value)) %>%
  ggplot(aes(fill=Source, y=value, x=new.names)) + 
  geom_bar(position="fill", stat="identity") + ylab("Relative influence") + xlab("Cluster") + theme_minimal() + 
  labs(title=expression(pi)) + theme(plot.title = element_text(size = 20)) + scale_fill_brewer(palette="Set1")




tlh_info_w <- as_tibble(lh_info_w)
colnames(tlh_info_w) <- paste0("Tree ", 1:10)
wlh_res<-pivot_longer(cbind(new.names,tlh_info_w), cols = 2:11)
wlh_res$Source <- c(rep("Likelihood", nrow(wlh_res)))

tp_info_w <- as_tibble(prior_info_w)
colnames(tp_info_w) <- paste0("Tree ", 1:10)
wp_res<-pivot_longer(cbind(new.names,tp_info_w), cols = 2:11)
wp_res$Source <- c(rep("Prior", nrow(wp_res)))

w_res <- rbind(wlh_res, wp_res)

w_res$new.names <- factor(w_res$new.names, levels = mixedsort(new.names))

# Stacked + percent. Average over the 10 trees
w_res <- w_res %>% group_by(new.names, name) %>% pivot_wider(names_from = "Source") %>% mutate(Likelihood = Likelihood/(Likelihood+Prior),
                                                                                      Prior = Prior/(Likelihood+Prior)) %>% 
  pivot_longer(cols = c("Likelihood", "Prior"), names_to = "Source")
  
  
p2 <- w_res %>% group_by(new.names, Source) %>% summarize(value = mean(value)) %>%
  ggplot(aes(fill=Source, y=value, x=new.names)) + 
  geom_bar(position="fill", stat="identity") + ylab("Relative influence") + xlab("Cluster") + theme_minimal() + 
  labs(title="w") + scale_fill_brewer(palette="Set1")


p1
p2

# Print to pdf
pdf(file = paste0("../Figures/", model, "/KL_divergence.pdf"), height = 8.3, width= 11.7, paper="a4r")
grid.arrange(p1, p2, nrow=2)
dev.off()
grid.arrange(p1, p2, nrow=2)



# Some likelihood plots
#like_Y[i, j] <- llOPTIM(pars2, si)

i <- seq(0,1,length=50)
j <- seq(0,1,length=50)
grid <- expand.grid(x=i,y=j)
pars <- unlist(init_)

# llOPTIM <- function(params, data)
grid <- ddply(grid,~x+y,mutate,loglik=llOPTIM(params=c(pars[1], pars[2], x, y),si))

grid <- subset(grid,is.finite(loglik))
ggplot(grid,aes(x=x,y=y,z=loglik,fill=loglik))+
  geom_tile()+geom_contour(binwidth=1)






