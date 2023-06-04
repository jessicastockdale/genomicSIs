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
mu <- 4.5
sigma <- 2.0 # si/gen time parameters
R0 <- 2
mutation_rate <- 1e-4
#

prop <- .5 #sampling proportion
low <- c(.01, .01, .01, .01); upp <- c(10, 10, .99, .99) # lower & upper bound
#dlt <- Inf; q <- 0.005
N <- 5 # how many outbreaks to generate
n <- 100 #of transmission trees to sample
#


# determine prior for pi, w
# + plot the prior distributions
Scenarios <- tibble(shape1 = c(12, 18, 6, 5 , 50, 9, 3), shape2 = c(12,12,12, 5,50, 3, 10))
# 1. 12, 12. Mean 0.5, SD .1 Baseline
# 2. 18, 12 Mean 0.6, SD 0.08. Inc mean
# 3. 6, 12 Mean 0.33, SD .1 Dec mean
# 4. 5, 5 Mean 0.5 SD 0.15. Inc SD
# 5. 50, 50 Mean 0.5, SD 0.05. Dec SD
Scen_names <- c("Baseline", "Increased mean", "Decreased mean", "Increased standard deviation", "Decreased standard deviation", "Further increased mean", "Further decreased mean")
# 1. larger mean, 2. smaller mean, 3. larger sd, 4. smaller sd, 5. baseline

# Plot the scenarios
x <- seq(0,1,length.out=1000)
scurves <- tibble(x = x, Baseline = dbeta(x,12,12), "Increased mean" = dbeta(x,18,12), "Decreased mean" = dbeta(x,6,12),
                  "Increased standard deviation" = dbeta(x,5,5), "Decreased standard deviation" = dbeta(x,50,50), "Further increased mean" = dbeta(x,9,3), "Further decreased mean" = dbeta(x,3,10))


mycols <- colorRampPalette(brewer.pal(8, "Spectral"))(7)

scurves %>% pivot_longer(!x, names_to = "Prior scenario") %>% 
  mutate(`Prior scenario` = factor(`Prior scenario`, c("Baseline", "Increased mean","Further increased mean", "Decreased mean","Further decreased mean", "Increased standard deviation","Decreased standard deviation"))) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_line(aes(color = `Prior scenario`), lwd=0.7) + theme_classic() + theme(legend.position="top", text = element_text(size=14)) +
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
  xlab(expression('Prior probability (' ~pi*' or '~w*')')) + ylab("Density") + scale_colour_manual(values = mycols)
ggsave("../simulation_study/prior_simscenarios.pdf", width = 9, height = 5.43, units = "in")
##

data.sum <- NULL

start <- Sys.time()

set.seed(57382)
for(j in 1:nrow(Scenarios)){
  for(i in 1:N){
    
    prior.pi <- unlist(Scenarios[j,])
    prior.w <- unlist(Scenarios[j,])
    
    init <- c(mu, sigma, prop, prop) # initial values of params
    
    # generate random time series and contact tracing
    outbreak <- simOutbreak2(nsus, c(mu, sigma), t, R0, mutation_rate, nsus/10, F, F)
    outbreak$epidata <- wiwAdj(outbreak$epidata, prop)
    #full <- createTransCloud(outbreak, max_genetic_dist = Inf)
    #eps <- quantile(full$distance, q)
    estim.cloud <- createTransCloud(outbreak, max_genetic_dist = mu*2/10000, onset_interval_diff = c(0, 35))
    
    # estimate the si from sampled trees
    trees <- sampleTransTree(estim.cloud, n)
    saveRDS(trees, file = paste0("../simulation_study/sampled_trees/",j, "number", i, ".rds"))
    estim.map <- siEstimate(trees, init, low, upp, .95, T, estimation_method = "map", config_setup = list(prior.pi = prior.pi, prior.w = prior.w))
    #estim.map <- siEstimate(trees, init, low, upp, .95, "no censoring", F, list(prior.pi = prior.pi))
    #
    
    par <- as.list(estim.map$estim$par)
    # Add CIs:
    CIpar <- estim.map$estim$bound
    par <- append(par, as.list(unmatrix(CIpar)))
    
    xx <- data.frame(par) %>% mutate(eps = mu*2/10000, Scen_name = Scen_names[j], pi.a = prior.pi[1], pi.b = prior.pi[2], w.a = prior.w[1], w.b = prior.w[2])
    data.sum <- bind_rows(data.sum, xx)
  }
}

write.csv(data.sum, "../simulation_study/sim_prior_res.csv")




