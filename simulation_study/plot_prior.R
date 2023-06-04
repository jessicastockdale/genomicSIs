library(ggplot2)
library(dplyr)
library(patchwork)
library(foreach)
library(readxl)
library(stringr)
library(tidyr)
library(gtools)

sumdat <- read.csv("../simulation_study/sim_prior_res.csv")

sumdat$X <- 1:nrow(sumdat)
sumdat$Scen_name <- factor(sumdat$Scen_name, levels = c(levels=c("Baseline", "Increased mean","Further increased mean", "Decreased mean","Further decreased mean", "Increased standard deviation","Decreased standard deviation")))

mycols <- colorRampPalette(brewer.pal(8, "Spectral"))(7)

# plot panel

pl.mu <- ggplot() + labs(x = "Prior scenario", y = expression(mu)) + theme_light() +
  #geom_jitter(aes(x = prop, y = mu), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = Scen_name, ymin = mu.lower, ymax = mu.upper, group=X, col = Scen_name), sumdat, width = .4, position=position_dodge(0.9)) +
  geom_point(aes(x = Scen_name, y = mu, group=X, col = Scen_name), sumdat, position=position_dodge(0.9)) +
  geom_hline(aes(yintercept = 4.5), lty = "dashed")+ theme(legend.position = "none") + scale_colour_manual(values = mycols) +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))

pl.sg <- ggplot() + labs(x = "Prior scenario", y = expression(sigma)) + theme_light() +
  #geom_jitter(aes(x = prop, y = sigma), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = Scen_name, ymin = sigma.lower, ymax = sigma.upper, group=X, col = Scen_name), sumdat,  width = .2, position=position_dodge(0.9)) +
  geom_point(aes(x = Scen_name, y = sigma, group=X, col = Scen_name), sumdat,  position=position_dodge(0.9)) + 
  geom_hline(aes(yintercept = 1), lty = "dashed")+theme(legend.position = "none") + scale_colour_manual(values = mycols)+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))

pl.pi <- ggplot() + labs(x = "Prior scenario", y = expression(pi)) + theme_light() +
  #geom_jitter(aes(x = prop, y = pi), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = Scen_name, ymin = pi.lower, ymax = pi.upper, group=X, col = Scen_name), sumdat,  width = .2, position=position_dodge(0.9)) +
  geom_point(aes(x = Scen_name, y = pi, group=X, col = Scen_name), sumdat,  position=position_dodge(0.9)) + theme(legend.position = "none") + scale_colour_manual(values = mycols)+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))
#geom_hline(aes(yintercept = 4.5), lty = "dashed")

pl.w <- ggplot() + labs(x = "Prior scenario", y = "w") + theme_light() +
  #geom_jitter(aes(x = prop, y = w), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = Scen_name, ymin = w.lower, ymax = w.upper, group=X, col = Scen_name), sumdat,  width = .2, position=position_dodge(0.9)) +
  geom_point(aes(x = Scen_name, y = w, group=X, col = Scen_name), sumdat, position=position_dodge(0.9)) + theme(legend.position = "none") + scale_colour_manual(values = mycols)+
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15))
#geom_hline(aes(yintercept = 4.5), lty = "dashed")

pl.mu/pl.sg/pl.pi/pl.w
ggsave("../simulation_study/prior_simests.pdf", width = 8, height = 7, units = "in")








# Do the KL analysis

library(LaplacesDemon)

# Structure for results
# 5 outbreaks in 5 scenarios (not 'further'), sample 5 of the trees
lh_info <- matrix(NA, 5*5, 5)
prior_info <- matrix(NA, 5*5, 5)


# Required functions
findpost_pi <- function(k, pars, priors){
  pars2 <- c(pars[1], pars[2], k, pars[4])
  # Post
  post_X <- -logll(pars2, si, .cop = T, mth = "map", ctr = list(prior.pi = priors, prior.w = priors))
}
findlike_pi <- function(k, pars, priors){
  pars2 <- c(pars[1], pars[2], k, pars[4])
  # Lh
  like_Y <- -logll(pars2, si, .cop = T, mth = "mle", ctr = list(prior.pi = priors, prior.w = priors))
  
}

findpost_w <- function(k, pars, priors){
  pars2 <- c(pars[1], pars[2], pars[3], k)
  # Post
  post_X <- -logll(pars2, si, .cop = T, mth = "map", ctr = list(prior.pi = priors, prior.w = priors))
}
findlike_w <- function(k, pars, priors){
  pars2 <- c(pars[1], pars[2], pars[3], k)
  # Lh
  like_Y <- -logll(pars2, si, .cop = T, mth = "mle", ctr = list(prior.pi = priors, prior.w = priors))
}


### Run it: PI PRIOR
for (i in 1:5){ # outbreak
  for (k in 1:5){ # scenario
    
    print(paste0("Currently running ", i, " in scenario ", k))
    
    # The beta prior: Beta( , )
    the_vec <- seq(0,1,0.001)
    prior_Y <- dbeta(the_vec, sumdat$pi.a[(i-1)*5+k], sumdat$pi.b[(i-1)*5+k])
    
    # Load sampled trees
    trees_sam <- readRDS(file = paste0("../simulation_study/sampled_trees/",k, "number", i, ".rds"))
    
    count <- 1
    
    # Sample 5 of the trees to compute
    for (j in sample(0:(length(trees_sam)-1))){
      
      if (count<6){
        
        print(paste0("Currently running tree", j))
        
        # Calculate posterior and likelihood densities
        si=trees_sam[[j+1]]
        si <- si[!is.na(si)]
        pars = init <- c(4.5, 1, 0.5, 0.5) # Will change pi/w to calculate the KL. Set remaining params to their maps
        pars[1] = sumdat$mu[(i-1)*5+k]
        pars[2] = sumdat$sigma[(i-1)*5+k]
        pars[4] = sumdat$w[(i-1)*5+k]
        
        post_X <- sapply(the_vec, findpost_pi, pars = pars, priors = c(sumdat$pi.a[(i-1)*5+k], sumdat$pi.b[(i-1)*5+k]))
        like_Y <- sapply(the_vec, findlike_pi, pars = pars, priors = c(sumdat$pi.a[(i-1)*5+k], sumdat$pi.b[(i-1)*5+k]))
        
        # Normalize
        post_X <- exp(post_X)/sum(exp(post_X))
        like_Y <- exp(like_Y)/sum(exp(like_Y))
        
        # Calculate KL divergence
        postprior <- KLD(post_X, prior_Y)
        postlike <- KLD(post_X, like_Y)
        lh_info[(i-1)*5+k,count] <- postprior$sum.KLD.px.py
        prior_info[(i-1)*5+k,count] <- postlike$sum.KLD.px.py
        print(paste0("lh info ", lh_info[(i-1)*5+k,count])) 
        print(paste0("prior info ", prior_info[(i-1)*5+k,count])) 
        
        
        count <- count + 1
        
        
        
      }
      
      plot(the_vec, prior_Y)
      plot(the_vec, like_Y)
      plot(the_vec, post_X)
      
    }
    
  }
}

lh_info_pi <- lh_info
prior_info_pi <- prior_info


# Plot the results (% stacked bars)

new.names = sumdat$Scen_name[1:25] # 5 outbreaks in 5 scenarios

tlh_info_pi <- as_tibble(lh_info_pi)
colnames(tlh_info_pi) <- paste0("Tree ", 1:5)
pilh_res<-pivot_longer(cbind(new.names,tlh_info_pi), cols = 2:6)
pilh_res$Source <- c(rep("Likelihood", nrow(pilh_res)))

tp_info_pi <- as_tibble(prior_info_pi)
colnames(tp_info_pi) <- paste0("Tree ", 1:5)
pip_res<-pivot_longer(cbind(new.names,tp_info_pi), cols = 2:6)
pip_res$Source <- c(rep("Prior", nrow(pip_res)))

pi_res <- rbind(pilh_res, pip_res)

pi_res$new.names <- factor(pi_res$new.names)

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








