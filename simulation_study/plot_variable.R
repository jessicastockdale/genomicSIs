library(ggplot2)
library(dplyr)
library(patchwork)
library(foreach)
library(readxl)
library(RColorBrewer)

sumdat <- read.csv("../simulation_study/sim_variable_res_set.csv")

sumdat$mu_true <- as.character(sumdat$mu_true)
sumdat$X <- 1:nrow(sumdat)

mutib <- tibble(mu = c(1.2, 2, 3, 4, 5), x = c(1,2,3,4,5))


# plot panel

mycols <- colorRampPalette(brewer.pal(8, "Spectral"))(6)

pl.mu <- ggplot() + labs(x = "True serial interval mean", y = expression(mu)) + theme_light() +
  #geom_jitter(aes(x = prop, y = mu), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = mu_true, ymin = mu.lower, ymax = mu.upper, group=X, col = mu_true), sumdat, width = .4, position=position_dodge(0.9)) +
  geom_point(aes(x = mu_true, y = mu, group=X, col = mu_true), sumdat, position=position_dodge(0.9)) + geom_point(aes(x = x, y=mu), mutib, size = 2) +
  #geom_hline(aes(yintercept = 4.5), lty = "dashed")+ 
  theme(legend.position = "none") + scale_colour_manual(values = mycols)
pl.mu

pl.sg <- ggplot() + labs(x = "True serial interval mean", y = expression(sigma)) + theme_light() +
  #geom_jitter(aes(x = prop, y = sigma), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = mu_true, ymin = sigma.lower, ymax = sigma.upper, group=X, col = mu_true), sumdat,  width = .2, position=position_dodge(0.9)) +
  geom_point(aes(x = mu_true, y = sigma, group=X, col = mu_true), sumdat,  position=position_dodge(0.9)) + 
  geom_hline(aes(yintercept = 1), lty = "dashed") + theme(legend.position = "none")  + scale_colour_manual(values = mycols)

pl.pi <- ggplot() + labs(x = "True serial interval mean", y = expression(pi)) + theme_light() +
  #geom_jitter(aes(x = prop, y = pi), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = mu_true, ymin = pi.lower, ymax = pi.upper, group=X, col = mu_true), sumdat,  width = .2, position=position_dodge(0.9)) +
  geom_point(aes(x = mu_true, y = pi, group=X, col = mu_true), sumdat,  position=position_dodge(0.9)) + theme(legend.position = "none")  + scale_colour_manual(values = mycols)
#geom_hline(aes(yintercept = 4.5), lty = "dashed")

pl.w <- ggplot() + labs(x = "True serial interval mean", y = "w") + theme_light() +
  #geom_jitter(aes(x = prop, y = w), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = mu_true, ymin = w.lower, ymax = w.upper, group=X, col = mu_true), sumdat,  width = .2, position=position_dodge(0.9)) +
  geom_point(aes(x = mu_true, y = w, group=X, col = mu_true), sumdat, position=position_dodge(0.9)) + theme(legend.position = "none") + scale_colour_manual(values = mycols)
#geom_hline(aes(yintercept = 4.5), lty = "dashed")

(pl.mu+pl.sg)/(pl.pi+pl.w)

ggsave("../simulation_study/mu_ests.pdf", width = 9, height = 5.43, units = "in")


# Relative error in mu and sd
sumdat$re_mu <- abs(sumdat$mu - as.numeric(sumdat$mu_true))/as.numeric(sumdat$mu_true)
sumdat$re_sigma <- abs(sumdat$sigma - 1)/1

# Stat test of significant difference


# do we think there are different means under each scenario (where they really do have different means)

mu = sumdat$mu
mu_se = (sumdat$mu.upper - sumdat$mu)/1.959964
Mu = sum((mu/mu_se^2)/sum(1/mu_se^2))
Stat=sum((mu - Mu)^2/mu_se)
Stat
# crit value chisq(24, 0.05) = 36.415
# Strong signal of differences amongst the 25 

mu = sumdat$mu[1:5]
mu_se = ((sumdat$mu.upper - sumdat$mu)/1.959964)[1:5]
Mu = sum((mu/mu_se^2)/sum(1/mu_se^2))
Stat=sum((mu - Mu)^2/mu_se)
Stat
# crit value chisq(4, 0.05) = 9.488
# no signif diff WITHIN 1st scenario - expected

mu = sumdat$mu[21:25]
mu_se = ((sumdat$mu.upper - sumdat$mu)/1.959964)[21:25]
Mu = sum((mu/mu_se^2)/sum(1/mu_se^2))
Stat=sum((mu - Mu)^2/mu_se)
Stat
# crit value chisq(4, 0.05) = 9.488
# same in the others


# Now look grouped by scenario
mu = sumdat$mu
mu_se = (sumdat$mu.upper - sumdat$mu)/1.959964
mu_typ <- 1
mu_se_typ <- 1
for (i in 1:5){
  mu_typ[i] = sum((mu[((i-1)*5+1):(i*5)]/mu_se[((i-1)*5+1):(i*5)]^2)/sum(1/mu_se[((i-1)*5+1):(i*5)]^2))
  mu_se_typ[i] = (sumdat$mu.upper[((i-1)*5+1):(i*5)] - sumdat$mu[((i-1)*5+1):(i*5)])/1.959964
}
Mu = sum((mu_typ/mu_se_typ^2)/sum(1/mu_se_typ^2))
Stat=sum((mu_typ - Mu)^2/mu_se_typ)
Stat
# crit value chisq(4, 0.05) = 9.488


mu = sumdat$mu
mu_se = (sumdat$mu.upper - sumdat$mu)/1.959964
mu_typ <- 1
mu_se_typ <- 1
for (i in 2:5){
  mu_typ[i] = sum((mu[((i-1)*5+1):(i*5)]/mu_se[((i-1)*5+1):(i*5)]^2)/sum(1/mu_se[((i-1)*5+1):(i*5)]^2))
  mu_se_typ[i] = (sumdat$mu.upper[((i-1)*5+1):(i*5)] - sumdat$mu[((i-1)*5+1):(i*5)])/1.959964
}
Mu = sum((mu_typ/mu_se_typ^2)/sum(1/mu_se_typ^2))
Stat=sum((mu_typ - Mu)^2/mu_se_typ)
Stat
# crit value chisq(3, 0.05) = 7.815



# w and pi
pi = sumdat$pi
pi_se = (sumdat$pi.upper - sumdat$pi)/1.959964
Mu = sum((pi/pi_se^2)/sum(1/pi_se^2))
Stat=sum((pi - Mu)^2/pi_se)
Stat

w = sumdat$w
w_se = (sumdat$w.upper - sumdat$w)/1.959964
Mu = sum((w/w_se^2)/sum(1/w_se^2))
Stat=sum((w - Mu)^2/w_se)
Stat
# crit value chisq(24, 0.05) = 36.415

mu = sumdat$pi
mu_se = (sumdat$pi.upper - sumdat$pi)/1.959964
mu_typ <- 1
mu_se_typ <- 1
for (i in 1:5){
  mu_typ[i] = sum((mu[((i-1)*5+1):(i*5)]/mu_se[((i-1)*5+1):(i*5)]^2)/sum(1/mu_se[((i-1)*5+1):(i*5)]^2))
  mu_se_typ[i] = (sumdat$mu.upper[((i-1)*5+1):(i*5)] - sumdat$mu[((i-1)*5+1):(i*5)])/1.959964
}
Mu = sum((mu_typ/mu_se_typ^2)/sum(1/mu_se_typ^2))
Stat=sum((mu_typ - Mu)^2/mu_se_typ)
Stat
# crit value chisq(4, 0.05) = 9.488

mu = sumdat$w
mu_se = (sumdat$w.upper - sumdat$w)/1.959964
mu_typ <- 1
mu_se_typ <- 1
for (i in 1:5){
  mu_typ[i] = sum((mu[((i-1)*5+1):(i*5)]/mu_se[((i-1)*5+1):(i*5)]^2)/sum(1/mu_se[((i-1)*5+1):(i*5)]^2))
  mu_se_typ[i] = (sumdat$mu.upper[((i-1)*5+1):(i*5)] - sumdat$mu[((i-1)*5+1):(i*5)])/1.959964
}
Mu = sum((mu_typ/mu_se_typ^2)/sum(1/mu_se_typ^2))
Stat=sum((mu_typ - Mu)^2/mu_se_typ)
Stat
# crit value chisq(4, 0.05) = 9.488