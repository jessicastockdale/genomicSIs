### Export results table for publication version

model <- if_else(coprim.transm==FALSE, "nc", "cop")
if (pi.model=="prior"){model <- paste0(model, "_priorpi")}
if (w.model=="prior"){model <- paste0(model, "w")}
model <- paste0(model, "_w", which.wave)

res <- readRDS(file = paste0(model, "_resultsobject.rda"))

res_tab <- tibble(cluster = sub(".*_", "", names), mu = 
  unlist(lapply(res$results, function(x){paste0(round(x$estimates$par[1], digits=2), " (", round(x$estimates$confint[1,1], digits=2) , ", ", round(x$estimates$confint[1,2], digits=2), ")")})),
  sigma = 
    unlist(lapply(res$results, function(x){paste0(round(x$estimates$par[2], digits=2), " (", round(x$estimates$confint[2,1], digits=2) , ", ", round(x$estimates$confint[2,2], digits=2), ")")})),
  pi = 
    unlist(lapply(res$results, function(x){paste0(round(x$estimates$par[3], digits=2), " (", round(x$estimates$confint[3,1], digits=2) , ", ", round(x$estimates$confint[3,2], digits=2), ")")})),
  w = 
    unlist(lapply(res$results, function(x){paste0(round(x$estimates$par[4], digits=2), " (", round(x$estimates$confint[4,1], digits=2) , ", ", round(x$estimates$confint[4,2], digits=2), ")")})),
)

# Sort rows numerically
res_tab <- res_tab[match(mixedsort(res_tab$cluster), res_tab$cluster, ),]

if (pool.trees){
  res_tab <- rbind(res_tab, c("Total",  paste0(round(res$results_pooled$estimates$par[1], digits=2), " (", round(res$results_pooled$estimates$confint[1,1], digits=2) , ", ", round(res$results_pooled$estimates$confint[1,2], digits=2), ")"), 
                              paste0(round(res$results_pooled$estimates$par[2], digits=2), " (", round(res$results_pooled$estimates$confint[2,1], digits=2) , ", ", round(res$results_pooled$estimates$confint[2,2], digits=2), ")"),
                              paste0(round(res$results_pooled$estimates$par[3], digits=2), " (", round(res$results_pooled$estimates$confint[3,1], digits=2) , ", ", round(res$results_pooled$estimates$confint[3,2], digits=2), ")"),
                              paste0(round(res$results_pooled$estimates$par[4], digits=2), " (", round(res$results_pooled$estimates$confint[4,1], digits=2) , ", ", round(res$results_pooled$estimates$confint[4,2], digits=2), ")")))
}



write.csv(res_tab, paste0("../Figures/", model, "_restab.csv"), row.names=FALSE)


