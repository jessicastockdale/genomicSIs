require(doParallel)
require(doSNOW)
require(tcltk)
require(dplyr)
#source("coputilities.R")

SIestim <- function(trees, si.config=list(par0, lower, upper, ci), method="OPTIM", nameSI="onset_diff", cores=1, howmany = length(trees), pi_prior=NULL, w_prior=NULL, progress.bars=TRUE){
  start <- Sys.time()
  singleSI <- function(treedata,pars,low,upp,ci,method, nm=c(nameSI)){
    si <- treedata[,nm[1]]
    si <- si[!is.na(si)]
    if(method=="OPTIM") res <- OPTIMestim(si=si, par0=pars, lower=low, upper=upp, ci=ci, pi_prior=pi_prior, w_prior=w_prior)
    else stop("Invalid method input. method=c('OPTIM').")
    condpar <- res$par
    condse <- res$se
    cov_ms <- res$varcov[2,1]
    out <- data.frame(mu=condpar[1], sigma=condpar[2], pi=condpar[3], w=condpar[4],
                      se.mu=condse[1], se.sigma=condse[2], cov.ms=cov_ms, se.pi=condse[3], se.w=condse[4],
                      LL=res$logLik, msg=res$message, stringsAsFactors=F)
    return(out)
  }
  
  #perform parallelism
  cl <- makeCluster(cores, type="SOCK")
  registerDoSNOW(cl)
  if (progress.bars==TRUE){
    pb <- tkProgressBar("test progress bar", "Iteration", 0, length(trees), 0)
    progress <- function(n){
      info <- sprintf("%1.0f%% done", n/length(trees)*100)
      setTkProgressBar(pb, n, sprintf("test (%s)", info), info)
    }
    opts <- list(progress=progress)
  }else{opts <- list()}
  estim <- foreach(i=1:length(trees), .combine=rbind, .packages=c("doParallel"), 
                   .export=c("randomSI","dcgg2","dgdd","dfgd","llOPTIM","postOPTIM","registerDoParallel","OPTIMestim","reltol"), 
                   .options.snow=opts, .errorhandling="pass", .verbose=T) %dopar% {
                     tryCatch(singleSI(treedata=trees[[i]], pars=si.config$par0, low=si.config$lower, upp=si.config$upper, 
                                       ci=si.config$ci, method=method),
                              error=function(e){
                                msg <- conditionMessage(e)
                                out <- data.frame(mu=NA, sigma=NA, pi=NA, w=NA, se.mu=NA, se.sigma=NA, se.pi=NA, se.w=NA,
                                                  LL=NA, msg=msg, stringsAsFactors=F)
                                return(out)
                              })
    
    }
  stopCluster(cl)
  
  #filter estimates that are exactly at the bound
  temp <- estim
  estim <- estim %>% filter(msg=="CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH" | (msg=="ESTIMATES ARE EXACTLY AT THE BOUND." & (pi == 1.0)))
  if(nrow(estim)<1) estim <- temp #return everything if even though the estimates are at bound or error.
  
  # Option to filter outlier trees - defined as outside  mean +- 3*sd.
  #temp2 <- estim
  # means
  #temp2 <- temp2[(which(temp2$alpha < (mean(estim$alpha) + 3*sd(estim$alpha)))),]
  #temp2 <- temp2[(which(temp2$alpha > (mean(estim$alpha) - 3*sd(estim$alpha)))),]
  #temp2 <- temp2[(which(temp2$beta < (mean(estim$beta) + 3*sd(estim$beta)))),]
  #temp2 <- temp2[(which(temp2$beta > (mean(estim$beta) - 3*sd(estim$beta)))),]
  #temp2 <- temp2[(which(temp2$pi < (mean(estim$pi) + 3*sd(estim$pi)))),]
  #temp2 <- temp2[(which(temp2$pi > (mean(estim$pi) - 3*sd(estim$pi)))),]
  # sds
  #temp2 <- temp2[(which(temp2$se.alpha < (mean(estim$se.alpha) + 3*sd(estim$se.alpha)))),]
  #temp2 <- temp2[(which(temp2$se.alpha > (mean(estim$se.alpha) - 3*sd(estim$se.alpha)))),]
  #temp2 <- temp2[(which(temp2$se.beta < (mean(estim$se.beta) + 3*sd(estim$se.beta)))),]
  #temp2 <- temp2[(which(temp2$se.beta > (mean(estim$se.beta) - 3*sd(estim$se.beta)))),]
  #temp2 <- temp2[(which(temp2$se.pi < (mean(estim$se.pi) + 3*sd(estim$se.pi)))),]
  #temp2 <- temp2[(which(temp2$se.pi > (mean(estim$se.pi) - 3*sd(estim$se.pi)))),]
  
  #estim <- temp2
  
  # Return only the first 'howmany' estimates
  estim <- estim[1:min(howmany,nrow(estim)),]
  
  #get total var and unconditional estimation
  uncond.mu <- mean(estim$mu)
  uncond.sigma <- mean(estim$sigma)
  uncond.pi <- mean(estim$pi)
  uncond.w <- mean(estim$w)
  
  totVar.mu <- var(estim$mu) + mean((estim$se.mu)^2)
  totVar.sigma <- var(estim$sigma) + mean((estim$se.sigma)^2)
  totVar.pi <- var(estim$pi) + mean((estim$se.pi)^2)
  totVar.w <- var(estim$w) + mean((estim$se.w)^2)
  
  #confidence interval
  z <- qnorm(si.config$ci+(1-si.config$ci)/2)
  matci <- matrix(nrow=4,ncol=2); colnames(matci) <- c("lower","upper"); rownames(matci) <- c("mu","sigma","pi","w")
  matci[,1] <- c(uncond.mu-z*sqrt(totVar.mu), uncond.sigma-z*sqrt(totVar.sigma), uncond.pi-z*sqrt(totVar.pi), uncond.w-z*sqrt(totVar.w))
  matci[,2] <- c(uncond.mu+z*sqrt(totVar.mu), uncond.sigma+z*sqrt(totVar.sigma), uncond.pi+z*sqrt(totVar.pi), uncond.w+z*sqrt(totVar.w))
  matci[,1][matci[,1]<0] <- 0
  matci[3:4,2][matci[3:4,2]>1] <- 1
  
  #outputs
  par <- c(uncond.mu, uncond.sigma, uncond.pi, uncond.w); names(par) <- c("mu","sigma","pi","w")
  sd <- c(sqrt(totVar.mu), sqrt(totVar.sigma), sqrt(totVar.pi), sqrt(totVar.w)); names(sd) <- c("mu","sigma","pi","w")
  out <- list(estimates=list(par=par, se=sd, confint=matci, runtime=Sys.time()-start), record=estim)
  gc()
  return(out)
}
