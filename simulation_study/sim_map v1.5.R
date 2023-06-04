# TITLE: SIMULATION USING MAP

# ABOUT:
# Have the simulation with various sample proportion and 
# use MAP.

# REQUIREMENTS
# @import: - estimateSI_v1.3.R in <Serial Interval> project
source("../simulation_study/estimateSI_v1.3.R")

#
#---------------------- %%%%%%%%%%%%% ----------------------->
#
library(ggplot2)
library(patchwork)
library(beeswarm)
#
# <-------- Outbreak parameters ---------------------------->
#
nsus <- 1000 # initial susceptibles
t <- 5000 # time duration of outbreak
mu <- 4.5; sigma <- 2 # si/gen time parameters
R0 <- 2
mutation_rate <- 1e-4
#

prop <- .4 #sampling proportion
low <- c(1, .1, .01, .01); upp <- c(10, 5, .99, .99) # lower & upper bound
init <- c(mu, sigma, prop, prop) # initial values of params
dlt <- Inf; q <- 0.005
N <- 1 # how many outbreaks to generate
n <- 1000 #of transmission trees to sample
#


# determine prior for pi
if(prop-.2 > 0) x1 <- prop-.2 else x1 <- .001
if(prop+.2 < 1) x3 <- prop+.2 else x3 <- .999

prior.pi <- findBeta(list(p = .5, x = prop), list(p = .99, x = x3), list(p = .01, x = x1))
prior.w <- findBeta(list(p = .5, x = .5), list(p = .99, x = .7), list(p = .01, x = .3))
#prior.w <- prior.pi

data.sum <- NULL

start <- Sys.time()


for(i in 1:N){
  # generate random time series and contact tracing
  set.seed(1234567)
  outbreak <- simOutbreak2(nsus, c(mu, sigma), t, R0, mutation_rate, 100, F, F)
  outbreak$epidata <- wiwAdj(outbreak$epidata, prop)
  full <- createTransCloud(outbreak, max_genetic_dist = Inf)
  eps <- quantile(full$distance, q)
  estim.cloud <- createTransCloud(outbreak, max_genetic_dist = eps)
  
  # estimate the si from sampled trees
  trees <- sampleTransTree(estim.cloud, n)
  estim.map <- siEstimate(trees, init, low, upp, .95, T, estimation_method = "map", config_setup = list(prior.pi = prior.pi, prior.w = prior.w))
  #estim.map <- siEstimate(trees, init, low, upp, .95, "no censoring", F, list(prior.pi = prior.pi))
  #
  
  par <- as.list(estim.map$estim$par)
  xx <- data.frame(par) %>% mutate(eps = eps, prop = prop, pi.a = prior.pi[1], pi.b = prior.pi[2], w.a = prior.w[1], w.b = prior.w[2])
  data.sum <- bind_rows(data.sum, xx)
}

