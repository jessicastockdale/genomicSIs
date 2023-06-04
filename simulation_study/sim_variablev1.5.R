# TITLE: SIMULATION USING MAP - varying SI variability

# REQUIREMENTS
# @import: - estimateSI_v1.3.R in <Serial Interval> project
source("../simulation_study/estimateSI_v1.3.R")

#
#---------------------- %%%%%%%%%%%%% ----------------------->
#
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(beeswarm)
library(gdata)
#
# <-------- Outbreak parameters ---------------------------->
#
nsus <- 200 # initial susceptibles
t <- 5000 # time duration of outbreak
mu <- c(1.2,2,3,4,5,6) # 4.5
sigma <- 1 # c(1,2,3,4,5) # si/gen time parameters
R0 <- 2
mutation_rate <- 1e-4
#

# Plot the distributions
x <- seq(0, 15, length.out = 100)
a <- mu^2/sigma^2; b <- mu/sigma^2

scurves <- tibble(x = x, "1.2" = dgamma(x,a[1],b[1]), "2.0" = dgamma(x,a[2],b[2]), "3.0" = dgamma(x,a[3],b[3]),
                  "4.0" = dgamma(x,a[4],b[4]), "5.0" = dgamma(x,a[5],b[5]))


mycols <- colorRampPalette(brewer.pal(8, "Spectral"))(6)

scurves %>% pivot_longer(!x, names_to = "Serial interval mean") %>% 
  mutate(`Serial interval mean` = factor(`Serial interval mean`, levels=colnames(scurves)[-1])) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_line(aes(color = `Serial interval mean`), lwd=1) + theme_classic() + theme(legend.position="top", text = element_text(size=14)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  xlab(expression('Time (days)')) + ylab("Density") + scale_colour_manual(values = mycols)
ggsave("../simulation_study/mu_scenarios.pdf", width = 9, height = 5.43, units = "in")
##

prop <- .5 #sampling proportion
low <- c(.01, .01, .01, .01); upp <- c(10, 10, .99, .99) # lower & upper bound
#dlt <- Inf; q <- 0.005
N <- 5 # how many outbreaks to generate
n <- 100 #of transmission trees to sample
#


# determine prior for pi
#if(prop-.2 > 0) x1 <- prop-.2 else x1 <- .001
#if(prop+.2 < 1) x3 <- prop+.2 else x3 <- .999

#prior.pi <- findBeta(list(p = .5, x = prop), list(p = .99, x = x3), list(p = .01, x = x1))
#prior.w <- findBeta(list(p = .5, x = .5), list(p = .99, x = .7), list(p = .01, x = .3))
prior.pi <- c(12,11)
prior.w <- prior.pi

data.sum <- NULL

start <- Sys.time()

for(j in 1:length(mu)){
  set.seed(8937878)
  
for(i in 1:N){
  
  init <- c(mu[j], sigma, prop, prop) # initial values of params
  
  # generate random time series and contact tracing
  outbreak <- simOutbreak2(nsus, c(mu[j], sigma), t, R0, mutation_rate, nsus/10, F, F)
  outbreak$epidata <- wiwAdj(outbreak$epidata, prop)
  #full <- createTransCloud(outbreak, max_genetic_dist = Inf)
  #eps <- quantile(full$distance, q)
  estim.cloud <- createTransCloud(outbreak, max_genetic_dist = mu[j]*2/10000, onset_interval_diff = c(0, 35)) #eps)
  
  # estimate the si from sampled trees
  trees <- sampleTransTree(estim.cloud, n)
  estim.map <- siEstimate(trees, init, low, upp, .95, T, estimation_method = "map", config_setup = list(prior.pi = prior.pi, prior.w = prior.w))
  #estim.map <- siEstimate(trees, init, low, upp, .95, "no censoring", F, list(prior.pi = prior.pi))
  #
  
  par <- as.list(estim.map$estim$par)
  # Add CIs:
  CIpar <- estim.map$estim$bound
  par <- append(par, as.list(unmatrix(CIpar)))
 
  mu_vec <- mu
  xx <- data.frame(par) %>% mutate(eps = mu[j]*2/10000, mu_true = mu_vec[j], pi.a = prior.pi[1], pi.b = prior.pi[2], w.a = prior.w[1], w.b = prior.w[2])
  data.sum <- bind_rows(data.sum, xx)
}
}

write.csv(data.sum, "../simulation_study/sim_variable_res_set.csv")

