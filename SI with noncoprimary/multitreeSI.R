require(doParallel)
require(doSNOW)
require(tcltk)
require(dplyr)
#source("cgg_utilities.R")

multitreeSI <- function(trees, si.config=list(init=list(mu,sigma,pi=NULL), lower, upper, ci), nameSI, cores=1, howmany = length(trees), pi_prior=NULL, progress.bars=TRUE){
  
  getSI <- function(data, nm=c(nameSI)){
    si <- data[,nm[1]]
    res <- siCGG(si=si[!is.na(si)], init=si.config$init, lower=si.config$lower, upper=si.config$upper, ci=si.config$ci,
                 pi_prior=pi_prior)
    condEst <- res$estimates
    condSe <- sqrt(abs(diag(res$varcov)))
    cov_ms <- res$varcov[2,1]
    out <- data.frame(mu=condEst[1], sigma=condEst[2], pi=condEst[3], se.mu=condSe[1], se.sigma=condSe[2],cov.ms=cov_ms, se.pi=condSe[3], LL=res$logLik, msg=res$message, stringsAsFactors=F)
    return(out)
  }
  
  #perform paralelism
  start <- Sys.time()
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
  estim <- foreach(i=1:length(trees), .combine=rbind, .packages=c("doParallel"), .export=c("dcgg","llcgg","lpcgg","registerDoParallel","siCGG"), .options.snow=opts) %dopar% {getSI(trees[[i]])}
  stopCluster(cl)
  
  #filter estimates that are exactly at the bound (alpha, beta) or otherwise give error
  temp <- estim
  estim <- estim %>% filter(msg=="CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"| (msg=="ESTIMATES ARE EXACTLY AT THE BOUND." & (pi == si.config$lower[3] | pi == si.config$upper[3])))
  if(nrow(estim)<1) estim <- temp #return everything if all the estimates are at bound or error.
  
  # Option to filter outlier trees - defined as outside  mean +- 3*sd.
  #temp2 <- estim
  # means
  #temp2 <- temp2[(which(temp2$mu < (mean(estim$mu) + 3*sd(estim$mu)))),]
  #temp2 <- temp2[(which(temp2$mu > (mean(estim$mu) - 3*sd(estim$mu)))),]
  #temp2 <- temp2[(which(temp2$sigma < (mean(estim$sigma) + 3*sd(estim$sigma)))),]
  #temp2 <- temp2[(which(temp2$sigma > (mean(estim$sigma) - 3*sd(estim$sigma)))),]
  #temp2 <- temp2[(which(temp2$pi < (mean(estim$pi) + 3*sd(estim$pi)))),]
  #temp2 <- temp2[(which(temp2$pi > (mean(estim$pi) - 3*sd(estim$pi)))),]
  # sds
  #temp2 <- temp2[(which(temp2$se.mu < (mean(estim$se.mu) + 3*sd(estim$se.mu)))),]
  #temp2 <- temp2[(which(temp2$se.mu > (mean(estim$se.mu) - 3*sd(estim$se.mu)))),]
  #temp2 <- temp2[(which(temp2$se.sigma < (mean(estim$se.sigma) + 3*sd(estim$se.sigma)))),]
  #temp2 <- temp2[(which(temp2$se.sigma > (mean(estim$se.sigma) - 3*sd(estim$se.sigma)))),]
  #temp2 <- temp2[(which(temp2$se.pi < (mean(estim$se.pi) + 3*sd(estim$se.pi)))),]
  #temp2 <- temp2[(which(temp2$se.pi > (mean(estim$se.pi) - 3*sd(estim$se.pi)))),]
  #estim <- temp2
  
  # Return only the first 'howmany' estimates
  estim <- estim[1:min(howmany,nrow(estim)),]
  
  #get total var and unconditional estimation
  uncond.mu <- mean(estim$mu)
  totVar.mu <- var(estim$mu) + mean((estim$se.mu)^2)
  uncond.sigma <- mean(estim$sigma)
  totVar.sigma <- var(estim$sigma) + mean((estim$se.sigma)^2)
  uncond.pi <- mean(estim$pi)
  totVar.pi <- var(estim$pi) + mean((estim$se.pi)^2)
  
  #get total covariance:
  ##formula Cov(a,b) = E(Cov(a,b|tau)) + cov(E(a|tau),E(b|tau))
  totCov.ms <- mean(estim$cov.ms) + cov(estim$mu, estim$sigma)
  
  #confidence interval
  z <- qnorm(si.config$ci+(1-si.config$ci)/2)
  matci <- matrix(nrow=3,ncol=2); colnames(matci) <- c("lower","upper"); rownames(matci) <- c("mu","sigma","pi")
  matci[,1] <- c(uncond.mu-z*sqrt(totVar.mu), uncond.sigma-z*sqrt(totVar.sigma), uncond.pi-z*sqrt(totVar.pi))
  matci[,2] <- c(uncond.mu+z*sqrt(totVar.mu), uncond.sigma+z*sqrt(totVar.sigma), uncond.pi+z*sqrt(totVar.pi))
  whichnegative <- which(matci[,1]<0 & !is.na(matci[,1]))
  if(length(whichnegative)>0) matci[,1][whichnegative] <- 0
  if(!is.na(matci[3,2]) && matci[3,2]>1) matci[3,2] <- 1

  print(Sys.time()-start)
  
  #outputs
  par <- c(uncond.mu,uncond.sigma,uncond.pi); names(par) <- c("mu","sigma","pi")
  sd <- c(sqrt(totVar.mu), sqrt(totVar.sigma), sqrt(totVar.pi)); names(sd) <- c("mu","sigma","pi")
  out <- list(estimates=list(par=par, se=sd, confint=matci), record=estim)
  gc()
  return(out)
}
