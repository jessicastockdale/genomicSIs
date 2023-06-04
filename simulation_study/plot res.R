library(ggplot2)
library(dplyr)
library(patchwork)
library(foreach)
library(readxl)

sum <- read_xlsx("./simulation/sim_map/result sim_map/result_allprop.xlsx")
rec <- read_xlsx("./simulation/sim_map/result sim_map/result_allprop.xlsx", 2)

rec$prop <- as.character(rec$prop)
sum$prop <- as.character(sum$prop)


# plot panel

pl.mu <- ggplot() + labs(x = "", y = expression(mu)) + theme_light() +
  #geom_jitter(aes(x = prop, y = mu), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = prop, ymin = low_mu, ymax = upp_mu), sum, col = "coral", width = .2) +
  geom_point(aes(x = prop, y = mu), sum, col = "coral") +
  geom_hline(aes(yintercept = 4.5), lty = "dashed")

pl.sg <- ggplot() + labs(x = "", y = expression(sigma)) + theme_light() +
  #geom_jitter(aes(x = prop, y = sigma), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = prop, ymin = low_sigma, ymax = upp_sigma), sum, col = "coral", width = .2) +
  geom_point(aes(x = prop, y = sigma), sum, col = "coral") +
  geom_hline(aes(yintercept = 2), lty = "dashed")

pl.pi <- ggplot() + labs(x = "", y = expression(pi)) + theme_light() +
  #geom_jitter(aes(x = prop, y = pi), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = prop, ymin = low_pi, ymax = upp_pi), sum, col = "coral", width = .2) +
  geom_point(aes(x = prop, y = pi), sum, col = "coral") 
  #geom_hline(aes(yintercept = 4.5), lty = "dashed")

pl.w <- ggplot() + labs(x = "", y = "w") + theme_light() +
  #geom_jitter(aes(x = prop, y = w), rec, col = "light blue", alpha = .4, width = .2) +
  geom_errorbar(aes(x = prop, ymin = low_w, ymax = upp_w), sum, col = "coral", width = .2) +
  geom_point(aes(x = prop, y = w), sum, col = "coral") 
  #geom_hline(aes(yintercept = 4.5), lty = "dashed")

(pl.mu+pl.sg)/(pl.pi+pl.w)


# plot distribution

x <- seq(0, 15, length.out = 100)
a <- sum$mu^2/sum$sigma^2; b <- sum$mu/sum$sigma^2
dt <- foreach(i = 1:10, .combine = "rbind") %do% data.frame(x = x, fx = dgamma(x, a[i], b[i]), prop = sum$prop[i])

pl.dist <- ggplot() + labs(x = "", y = "density") + theme_light() +
  guides(col = guide_legend(title = "")) +
  geom_line(aes(x = x, y = fx, col = prop), dt) +
  geom_line(aes(x = x, y = fx), data.frame(x = x, fx = dgamma(x, 4.5^2/4, 4.5/4)), lwd = 1.25)






# perform linear regression between pi and w
fit <- lm(w~pi, data = sum)
#summary(fit)


# plot it

pl.piw <- ggplot() + labs(x = expression(pi), y = "w") + theme_light() +
  #geom_point(aes(x = pi, y = w, col = prop), rec, alpha = .4) +
  geom_point(aes(x = pi, y = w), sum, shape = 4) +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], lty = "dashed", alpha = .2) +
  annotate(geom = "text", x = .9, y = .43, label = "intercept = 0.292") +
  annotate(geom = "text", x = .9, y = .4, label = "slope = 0.593")


x <- seq(0, 1, length.out = 100)
dt.pi <- foreach(i = 1:10, .combine = "rbind") %do% data.frame(x = x, fx = dbeta(x, sum$pi.a[i], sum$pi.b[i]), prop = sum$prop[i])
dt.w <- foreach(i = 1:10, .combine = "rbind") %do% data.frame(x = x, fx = dbeta(x, sum$w.a[i], sum$w.b[i]), prop = sum$prop[i])

hist.pi <- ggplot(rec) + labs(x = "", y = expression(paste("density of ", pi))) + theme_light() +
  geom_histogram(aes(x = pi, y = ..density..), rec, bins = 50, fill = "white", col = "black") +
  geom_line(aes(x = x, y = fx), dt.pi) +
  facet_wrap(~prop)
hist.w <- ggplot(rec) + labs(x = "", y = paste("density of w")) + theme_light() +
  geom_histogram(aes(x = w, y = ..density..), rec, bins = 50, fill = "white", col = "black") +
  geom_line(aes(x = x, y = fx), dt.w) +
  facet_wrap(~prop)

