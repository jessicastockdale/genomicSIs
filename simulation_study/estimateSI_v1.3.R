# TITLE: MAIN SCRIPT TO ESTIMATE SERIAL INTERVAL

# <---------------- %%%%%%%%%%%% ----------------------------------------->
#
require(doSNOW)
require(tcltk)
require(dplyr)
require(foreach)
require(igraph)
require(adegenet)
require(ape)
require(LearnBayes)
require(dfoptim)
#


#
# <---------------- Utilities functions ------------------------------->
#
rcgg <- function(n,mu,sigma,pi){
  # use shape and rate
  a <- (mu/sigma)^2; b <- mu/sigma^2
  m <- rgeom(n,pi)
  x <- rgamma(n, (m+1)*a, b)
  return(x)
}
dcgg <- function(x, mu, sigma, pi){
  gt <- function(xi,mu,sigma,pi){
    # use shape and rate
    a <- (mu/sigma)^2; b <- mu/sigma^2
    
    #Find maximum m within tol=1e-10
    k <- 0; tol <- 1e-10
    y <- 1/gamma(a)
    v <- b^a * xi^a * (1-pi)
    while(y>tol && y<Inf){
      k <- k+1
      y <- tryCatch(v^k/gamma((k+1)*a), error=function(e) return(0), warning=function(w) return(0))
      if(is.nan(y)) y <- 0 #break the loop
      #NaN happens when gamma((m+1)*a)=Inf or when it's not defined
    }
    if(k>1) mmax <- k-1
    else mmax <- 0
    f <- sapply(0:mmax, function(i) dgamma(xi, (i+1)*a, b) * dgeom(i, pi))
    return(sum(f))
  }
  ft <- sapply(1:length(x), function(i) gt(xi=x[i], mu, sigma, pi))
  ft[is.infinite(ft)] <- 0 # it evaluate x=0 as Inf. 
  return(ft)
}
pcgg <- function(x, mu, sigma, pi){
  pt <- function(xi,mu,sigma,pi){
    int <- integrate(f=dcgg, lower=0, upper=xi, mu=mu, sigma=sigma, pi=pi)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], mu, sigma, pi))
  Ft[which(Ft>1)] <- 1
  return(Ft)
}
rfgd <- function(n,mu,sigma){
  # use shape and rate
  a <- (mu/sigma)^2; b <- mu/sigma^2
  
  x1 <- rgamma(n, a, b)
  x2 <- rgamma(n, a, b)
  return(abs(x1-x2))
}
dfgd <- function(x, mu, sigma){
  gt <- function(ti, mu, sigma){
    a <- (mu/sigma)^2; b <- mu/sigma^2
    if(ti>=0){
      integrand <- function(u, t) dgamma(u,a,b)*dgamma(u-t,a,b)
      ht <- tryCatch(integrate(integrand, ti, Inf, t = ti)$value,
                     error = function(e) return(0))
    } else{
      integrand <- function(u, t) dgamma(u,a,b)*dgamma(u+t,a,b)
      ht <- tryCatch(integrate(integrand, -ti, Inf, t = ti)$value,
                     error = function(e) return(0))
    }
    if(is.nan(ht)) ht <- 0 #happens because we may multiply 0 and Inf
    return(ht)
  }
  ft <- sapply(1:length(x), function(i) gt(ti=x[i], mu, sigma))
  return(2*ft)
}
pfgd <- function(x,mu,sigma){
  pt <- function(xi,mu,sigma){
    int <- integrate(f=dfgd, lower=0, upper=xi, mu=mu, sigma=sigma)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], mu, sigma))
  Ft[which(Ft>1)] <- 1
  return(Ft)
}
logll <- function(params, dt, .cop, mth, ctr){
  
  # if the si is only non-coprimary
  # only allowed method: trunc., mle, map
  if(!.cop){
    
    # setup
    # get error when computing gamma(0)
    if(params[1] == 0) params[1] <- 1e-300 
    if(params[2] == 0) params[2] <- 1e-300
    if(params[3] == 0) params[3] <- 1e-300
    if(params[3] > 1) params[3] <- 1
    dt <- dt[dt > 0] #omit si=0 when using non-coprimary
    
    # if the method is neither mle, trunc., map
    if(!(mth %in% c("mle","map","trunc."))) mth <- "mle"
    
    # if the method is mle
    if(mth == "mle"){
      # compute log-likelihood
      ft <- dcgg(dt, params[1], params[2], params[3])
      ft <- ft[!is.na(ft)] # avoid NaN
      ft[ft==0 | is.infinite(ft)] <- 1e-300 # avoid log(0). L-BFGS-B only allows finite obj. func.
      nll <- -sum(log(ft))
    } #end mle
    
    # if the method is truncation
    if(mth == "trunc."){
      # check if all controlled params are provided
      if(is.null(ctr$epsilon0) || ctr$epsilon0 <= 0) ctr$epsilon0 <- 1
      
      # truncate the data based on epsilon0
      dt <- dt[dt >= ctr$epsilon0]
      
      # compute the log-ll
      ptrunc <- 1-pcgg(ctr$epsilon0, params[1], params[2], params[3])
      ft <- dcgg(dt, params[1], params[2], params[3])
      ft <- ft[!is.na(ft)]
      ft[ft==0 | is.infinite(ft)] <- 1e-300 # avoid log(0). L-BFGS-B only allows finite obj. func.
      nll <- -sum(log(ft)) + length(dt)*log(ptrunc)
    } #end truncation method
    
    # if the method is map
    if(mth == "map"){
      # check if all controlled params are provided
      if(is.null(ctr$prior.pi) || min(ctr$prior.pi) <= 0) ctr$prior.pi <- c(5,3) 
      
      # compute the log-ll
      ft <- dcgg(dt, params[1], params[2], params[3])
      fprior <- dbeta(params[3], ctr$prior.pi[1], ctr$prior.pi[2])
      ft <- ft[!is.na(ft)]
      ft[ft==0 | is.infinite(ft)] <- 1e-300 
      if(fprior == 0) fprior <- 1e-300
      nll <- -sum(log(ft)) - length(dt)*log(fprior)
    }
  } #end non-coprimary
  
  
  # if the si includes coprimary
  # only allowed method: mle, map
  if(.cop){
    
    #setup
    #get error when computing gamma(0)
    if(params[1] == 0) params[1] <- 1e-300 
    if(params[2] == 0) params[2] <- 1e-300
    if(params[3] == 0) params[3] <- 1e-300
    if(params[3] > 1) params[3] <- 1
    if(params[4] < 0) params[4] <- 0
    if(params[4] > 1) params[4] <- 1
    
    # if the method is neither mle nor map
    if(!(mth %in% c("mle","map"))) mth <- "mle" 
    
    # if method == mle
    if(mth == "mle"){
      f1 <- dcgg(dt, params[1], params[2], params[3])
      f2 <- dfgd(dt, params[1], params[2])
      f1 <- f1[!is.na(f1)]
      f2 <- f2[!is.na(f2)]
      f1[f1==0 | is.infinite(f1)] <- 1e-300 #avoid log(0)
      f2[f2==0 | is.infinite(f2)] <- 1e-300 #avoid log(0)
      #f1[is.infinite(f1)] <- 1e300
      #f2[is.infinite(f2)] <- 1e300
      nll <- -sum(log(params[4]*f1 + (1-params[4])*f2))
    }
    
    # if method == "map"
    if(mth == "map"){
      # check if all controlled params are provided
      if(is.null(ctr$prior.pi)) fp <- 1
      else fp <- dbeta(params[3], ctr$prior.pi[1], ctr$prior.pi[2])
      
      if(is.null(ctr$prior.w)) fw <- 1
      else fw <- dbeta(params[4], ctr$prior.w[1], ctr$prior.w[1])
      
      # compute the log-ll
      f1 <- dcgg(dt, params[1], params[2], params[3])
      f2 <- dfgd(dt, params[1], params[2])
      f1 <- f1[!is.na(f1)]
      f2 <- f2[!is.na(f2)]
      f1[f1==0 | is.infinite(f1)] <- 1e-300 
      f2[f2==0 | is.infinite(f2)] <- 1e-300
      #f1[is.infinite(f1)] <- 1e300
      #f2[is.infinite(f2)] <- 1e300
      if(fp == 0) fp <- 1e-300
      if(fw == 0) fw <- 1e-300
      nll <- -sum(log(params[4]*f1 + (1-params[4])*f2)) - log(fp) - log(fw)
    }
  } #end coprimary
  
  return(nll)
}
opt <- function(params0, dt, lower, upper, .cop, mth, ctr){
  
  # omit NA
  dt <- dt[!is.na(dt)]
  
  # check completeness of control
  if(!.cop){
    if(!(mth %in% c("mle","map","trunc."))) cat("method must be mle, map, or trunc. By default, pick mle.", "\n")
    if(mth == "trunc." && (is.null(ctr$epsilon0) || ctr$epsilon0 <= 0)) cat("epsilon0 must be >0. By default, pick epsilon0 = 1", "\n")
    if(mth == "map" && (is.null(ctr$prior.pi) || min(ctr$prior.pi) <= 0)) cat("prior.pi must be shape1>0 shape2>0. By default, pick (5, 3)", "\n")
  }
  if(.cop){
    if(!(mth %in% c("mle","map"))) cat("method must be mle or map. By default pick mle.", "\n")
    if(mth == "map" && (is.null(ctr$prior.pi) || min(ctr$prior.pi) <= 0)) cat("prior.pi must be shape1>0 shape2>0. By default, pick (5, 3)", "\n")
    if(mth == "map" && (is.null(ctr$prior.w) || min(ctr$prior.w) <= 0)) cat("prior.w must be shape1>0 shape2>0. By default, pick (5, 3)", "\n")
  }
  
  # estimate parameters
  #myestim <- optim(params0, logll, dt=dt, method="L-BFGS-B", ctr = ctr,
  #                 lower=lower, upper=upper, .cop = .cop, mth = mth)
  myestim <- dfoptim::nmkb(params0, logll, dt=dt, ctr = ctr, .cop = .cop, mth = mth,
                  lower=lower, upper=upper)
  
  params0 <- myestim$par
  mylogll <- myestim$value
  convergence <- myestim$convergence
  if(all(params0[1:2] > lower[1:2]) && all(params0[1:2] < upper[1:2])) msg <- myestim$message
  else{
    msg <- "ESTIMATES EXCEED BOUNDS."
    convergence <- -1
  }
  
  # estimate var-cov matrix using hessian
  hessian <- numDeriv::hessian(logll, params0, dt = dt, .cop = .cop, mth = mth, ctr = ctr)
  varcov <- solve(hessian)
  
  param <- rep(NA, 4); names(param) <- c("mu", "sigma", "pi", "w")
  se <- rep(NA, 4); names(se) <- c("mu", "sigma", "pi", "w")
  if(!.cop){
    param[1:3] <- myestim$par
    se[1:3] <- sqrt(abs(diag(varcov)))
  } 
  else{
    param[1:4] <- myestim$par
    se[1:4] <- sqrt(abs(diag(varcov)))
  } 
  
  out <- list(par=param, se=se, varcov=varcov, logll=mylogll, convergence = convergence, msg=msg)
  return(out)
}
findBeta <- function(quantile1, quantile2, quantile3){
  # find the quantiles specified by quantile1 and quantile2 and quantile3
  quantile1_p <- quantile1[[1]]; quantile1_q <- quantile1[[2]]
  quantile2_p <- quantile2[[1]]; quantile2_q <- quantile2[[2]]
  quantile3_p <- quantile3[[1]]; quantile3_q <- quantile3[[2]]
  
  # find the beta prior using quantile1 and quantile2
  priorA <- beta.select(quantile1,quantile2)
  priorA_a <- priorA[1]; priorA_b <- priorA[2]
  
  # find the beta prior using quantile1 and quantile3
  priorB <- beta.select(quantile1,quantile3)
  priorB_a <- priorB[1]; priorB_b <- priorB[2]
  
  # find the best possible beta prior
  diff_a <- abs(priorA_a - priorB_a); diff_b <- abs(priorB_b - priorB_b)
  step_a <- diff_a / 100; step_b <- diff_b / 100
  if (priorA_a < priorB_a) { start_a <- priorA_a; end_a <- priorB_a }
  else                     { start_a <- priorB_a; end_a <- priorA_a }
  if (priorA_b < priorB_b) { start_b <- priorA_b; end_b <- priorB_b }
  else                     { start_b <- priorB_b; end_b <- priorA_b }
  steps_a <- seq(from=start_a, to=end_a, length.out=1000)
  steps_b <- seq(from=start_b, to=end_b, length.out=1000)
  max_error <- 10000000000000000000
  best_a <- 0; best_b <- 0
  for (a in steps_a)
  {
    for (b in steps_b)
    {
      # priorC is beta(a,b)
      # find the quantile1_q, quantile2_q, quantile3_q quantiles of priorC:
      priorC_q1 <- qbeta(c(quantile1_p), a, b)
      priorC_q2 <- qbeta(c(quantile2_p), a, b)
      priorC_q3 <- qbeta(c(quantile3_p), a, b)
      priorC_error <- abs(priorC_q1-quantile1_q) +
        abs(priorC_q2-quantile2_q) +
        abs(priorC_q3-quantile3_q)
      if (priorC_error < max_error)
      {
        max_error <- priorC_error; best_a <- a; best_b <- b
      }
    }
  }
  #print(paste("The best beta prior has a=",best_a,"b=",best_b))
  return(c(best_a, best_b))
}
siEstimate <- function(sampled_trees,
                       init_parameters,
                       lower_bound,
                       upper_bound,
                       confidence_interval,
                       .useCoprimary = FALSE,
                       .showProgressBar = TRUE,
                       estimation_method,
                       config_setup){
  start <- Sys.time()
  
  #estimate one tree
  tmpfn <- function(params0, dt, low, upp, .cop, mth, ctr){
    dt <- dt[!is.na(dt)]
    result <- opt(params0, dt, low, upp, .cop, mth, ctr)
    
    # get all outputs
    par <- result$par
    se <- result$se
    output <- data.frame(mu = par[1], sigma = par[2], pi = par[3], w = par[4],
                         se.mu = se[1], se.sigma = se[2], se.pi = se[3], se.w = se[4],
                         ll = result$logll, convergence = result$convergence, msg = result$msg, stringsAsFactors = FALSE)
    return(output)
  }
  
  # perform parallel
  if(is.null(config_setup$ncore)) ncore <- parallel::detectCores()
  cl <- makeCluster(ncore, type="SOCK")
  registerDoSNOW(cl)
  if(.showProgressBar){
    pb <- tkProgressBar(paste("SI estimation using ", estimation_method), "Progress...", 0, length(sampled_trees), 0)
    progress <- function(n){
      info <- sprintf("%1.0f%% done", n/length(sampled_trees)*100)
      setTkProgressBar(pb, n, paste("SI estimation using ", estimation_method), info)
    }
    opts <- list(progress=progress)
    result <- foreach(i = 1:length(sampled_trees),
                      .combine = rbind,
                      .packages = c("doParallel", "parallel"),
                      .export = c("makeCluster", "dcgg", "dfgd", "logll", "opt", "pcgg"),
                      .options.snow=opts, .errorhandling="pass") %dopar% {
                        tryCatch(tmpfn(init_parameters, sampled_trees[[i]], lower_bound, upper_bound, .useCoprimary, estimation_method, config_setup),
                                 error = function(e){
                                   msg <- conditionMessage(e)
                                   output <- data.frame(mu = NA, sigma = NA, pi = NA, w = NA,
                                                        se.mu = NA, se.sigma = NA, se.pi = NA, se.w = NA,
                                                        ll = NA, convergence = NA, msg = msg, stringsAsFactors = FALSE)
                                   return(output)
                                 })
                      }
    close(pb)
  } else{
    result <- foreach(i = 1:length(sampled_trees),
                      .combine = rbind,
                      .packages = c("doParallel", "parallel"),
                      .export = c("makeCluster", "dcgg", "dfgd", "logll", "opt", "pcgg"),
                      .errorhandling="pass") %dopar% {
                        tryCatch(tmpfn(init_parameters, sampled_trees[[i]], lower_bound, upper_bound, .useCoprimary, estimation_method, config_setup),
                                 error = function(e){
                                   msg <- conditionMessage(e)
                                   output <- data.frame(mu = NA, sigma = NA, pi = NA, w = NA,
                                                        se.mu = NA, se.sigma = NA, se.pi = NA, se.w = NA,
                                                        ll = NA, convergence = NA, msg = msg, stringsAsFactors = FALSE)
                                   return(output)
                                 })
                      }
    
  }
  stopCluster(cl)
  
  ## NOW CREATE THE OUTPUTS
  #filter estimates that are exactly at the bound
  temp <- result
  #resultBound <- result %>% filter(convergence == -1)
  #estim <- estim %>% filter(!(msg %in% c("ESTIMATES ARE EXACTLY AT THE BOUND.", "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH")))
  result <- result %>% 
    filter(convergence == 0)
  if(nrow(result) < 1) result <- temp #return everything if even though the estimates are at bound or error.
  
  #get total var and unconditional estimation
  param <- rep(NA, 4); names(param) <- c("mu", "sigma", "pi", "w")
  se <- rep(NA, 4); names(se) <- c("mu", "sigma", "pi", "w")
  ci <- matrix(NA, ncol = 2, nrow = 4); rownames(ci) <- c("mu", "sigma", "pi", "w"); colnames(ci) <- c("lower", "upper")
  
  #confidence interval
  z <- qnorm(confidence_interval + (1 - confidence_interval)/2)
  
  sc <- nrow(result)/length(sampled_trees)
  
  if(nrow(result)==0) result <- temp
  if(!.useCoprimary){
    param[1:3] <- c(mean(result$mu), mean(result$sigma), mean(result$pi))
    se[1:3] <- sqrt(c(var(result$mu) + mean((result$se.mu)^2), var(result$sigma) + mean((result$se.sigma)^2), var(result$pi) + mean((result$se.pi)^2)))
    ci[1:3,1] <- c(param[1] - z*se[1], param[2] - z*se[2], param[3] - z*se[3])
    ci[1:3,2] <- c(param[1] + z*se[1], param[2] + z*se[2], param[3] + z*se[3])
    ci[1:3,1] <- ifelse(ci[1:3, 1] < 0, 0, ci[1:3,1])
    ci[3,2] <- ifelse(ci[3, 2] > 1, 1, ci[3,2])
  } else{
    param[1:4] <- c(mean(result$mu), mean(result$sigma), mean(result$pi), mean(result$w))
    se[1:4] <- sqrt(c(var(result$mu) + mean((result$se.mu)^2), var(result$sigma) + mean((result$se.sigma)^2), var(result$pi) + mean((result$se.pi)^2), var(result$w) + mean((result$se.w)^2)))
    ci[,1] <- c(param[1] - z*se[1], param[2] - z*se[2], param[3] - z*se[3], param[4] - z*se[4])
    ci[,2] <- c(param[1] + z*se[1], param[2] + z*se[2], param[3] + z*se[3], param[4] + z*se[4])
    ci[,1] <- ifelse(ci[, 1] < 0, 0, ci[,1])
    ci[3:4,2] <- ifelse(ci[3:4, 2] > 1, 1, ci[3:4,2])
  } 
  output <- list(estim = list(par = param, se = se, bound = ci, run.time = Sys.time()-start, success.run = sc), 
                 record = temp)
  return(output)
}
#

#
# <---------------- Simulate outbreak using outbreaker ------------------->
simOutbreak2 <- function(initial_sus = 200,
                         genint_params, # vector in mean, sigma
                         outbreak_duration = 100,
                         R0,
                         mutation_rate = 1e-4,
                         minimum_case = 10,
                         .plot = FALSE,
                         .export = FALSE){
  ##FROM OUTBREAKER
  theoutbreak <- function(R0, infec.curve, n.hosts=200, duration=50,
                          seq.length=1e4, mu.transi=1e-4, mu.transv=mu.transi/2,
                          rate.import.case=0.01, diverg.import=10, group.freq=1,
                          stop.once.cleared=TRUE){
    ## HANDLE ARGUMENTS ##
    ## handle group sizes
    if(any(group.freq<0)) stop("negative group frequencies provided")
    group.freq <- group.freq/sum(group.freq)
    K <- length(group.freq)
    ## host.group <- sample(1:K, size=n.hosts, prob=group.freq, replace=TRUE)
    R0 <- rep(R0, length=K) # recycle R0
    
    ## normalize gen.time
    infec.curve <- infec.curve/sum(infec.curve)
    infec.curve <- c(infec.curve, rep(0, duration)) # make sure dates go all the way
    t.clear <- which(diff(infec.curve<1e-10)==1) # time at which infection is cleared
    
    ## GENETIC FUNCTIONS ##
    NUCL <- as.DNAbin(c("a","t","c","g"))
    TRANSISET <- list('a'=as.DNAbin('g'), 'g'=as.DNAbin('a'), 'c'=as.DNAbin('t'), 't'=as.DNAbin('c'))
    TRANSVSET <- list('a'=as.DNAbin(c('c','t')), 'g'=as.DNAbin(c('c','t')), 'c'=as.DNAbin(c('a','g')), 't'=as.DNAbin(c('a','g')))
    
    
    ## AUXILIARY FUNCTIONS ##
    ## generate sequence from scratch
    seq.gen <- function(){
      ##res <- list(sample(NUCL, size=seq.length, replace=TRUE)) # DNAbin are no longer lists by default
      res <- sample(NUCL, size=seq.length, replace=TRUE)
      class(res) <- "DNAbin"
      return(res)
    }
    
    ## create substitutions for defined SNPs - no longer used
    substi <- function(snp){
      res <- sapply(1:length(snp), function(i) sample(setdiff(NUCL,snp[i]),1)) # ! sapply does not work on DNAbin vectors directly
      class(res) <- "DNAbin"
      return(res)
    }
    
    ## create transitions for defined SNPs
    transi <- function(snp){
      res <- unlist(TRANSISET[as.character(snp)])
      class(res) <- "DNAbin"
      return(res)
    }
    
    ## create transversions for defined SNPs
    transv <- function(snp){
      res <- sapply(TRANSVSET[as.character(snp)],sample,1)
      class(res) <- "DNAbin"
      return(res)
    }
    
    ## duplicate a sequence (including possible mutations)
    seq.dupli <- function(seq, T){ # T is the number of time units between ancestor and decendent
      ## transitions ##
      n.transi <- rbinom(n=1, size=seq.length*T, prob=mu.transi) # total number of transitions
      if(n.transi>0) {
        idx <- sample(1:seq.length, size=n.transi, replace=FALSE)
        seq[idx] <- transi(seq[idx])
      }
      
      ## transversions ##
      n.transv <- rbinom(n=1, size=seq.length*T, prob=mu.transv) # total number of transitions
      if(n.transv>0) {
        idx <- sample(1:seq.length, size=n.transv, replace=FALSE)
        seq[idx] <- transv(seq[idx])
      }
      return(seq)
    }
    
    ## define the group of 'n' hosts
    choose.group <- function(n){
      out <- sample(1:K, size=n, prob=group.freq, replace=TRUE)
      return(out)
    }
    
    
    ## MAIN FUNCTION ##
    ## initialize results ##
    dynam <- data.frame(nsus=integer(duration+1), ninf=integer(duration+1), nrec=integer(duration+1))
    rownames(dynam) <- 0:duration
    res <- list(n=1, dna=NULL, onset=NULL, id=NULL, ances=NULL, dynam=dynam)
    res$dynam$nsus[1] <- n.hosts-1
    res$dynam$ninf[1] <- 1
    res$onset[1] <- 0
    res$id <- 1 # id of infected individuals
    res$ances <- NA
    res$group <- choose.group(1)
    EVE <- seq.gen()
    res$dna <- matrix(seq.dupli(EVE, diverg.import),nrow=1)
    class(res$dna) <- "DNAbin"
    res$status <- c("I", rep("S", n.hosts-1)) # will be I, S, or R
    
    
    ## run outbreak ##
    for(t in 1:duration){
      ## DETERMINE NEW INTERNAL INFECTIONS ##
      ## individual force of infection - purely based on symptom onset
      indivForce <- infec.curve[t-res$onset+1]
      
      ## temporal (spatial) force of infection * R0
      indivForce <- indivForce * R0[res$group]
      
      ## global force of infection (R0 \sum_j I_t^j / N)
      N <- res$dynam$nrec[t] + res$dynam$ninf[t] + res$dynam$nsus[t] # this may change because of imports
      globForce <- sum(indivForce)/N
      
      ## stop if no ongoing infection in the population
      if(stop.once.cleared && (globForce < 1e-12)) break;
      
      ## compute proba of infection for each susceptible
      p <- 1-exp(-globForce)
      
      ## number of new infections
      nbNewInf <- rbinom(1, size=res$dynam$nsus[t], prob=p)
      
      
      ## HANDLE NEW INTERNAL INFECTIONS ##
      if(nbNewInf>0){
        ## dates of new infections ##
        res$onset <- c(res$onset, rep(t,nbNewInf))
        
        ## identify the infectors of the new cases ##
        newAnces <- sample(res$id, size=nbNewInf, replace=TRUE, prob=indivForce)
        res$ances <- c(res$ances,newAnces)
        
        ## find the groups of the new cases ##
        newGroup <- choose.group(nbNewInf)
        res$group <- c(res$group,newGroup)
        
        ## id of the new cases ##
        areSus <- which(res$status=="S") # IDs of susceptibles
        newId <- sample(areSus, size=nbNewInf, replace=FALSE)
        res$id <- c(res$id, newId)
        res$status[newId] <- "I"
        
        ## dna sequences of the new cases ##
        ## molecular clock / generation
        ## newSeq <- t(sapply(match(newAnces, res$id), function(i) seq.dupli(res$dna[i,], 1)))
        ## molecular clock / time unit
        newSeq <- t(sapply(match(newAnces, res$id), function(i) seq.dupli(res$dna[i,], t-res$onset[match(newAnces, res$id)])))
        res$dna <- rbind(res$dna, newSeq)
      }
      
      
      ## IMPORTED CASES ##
      ## number of imported cases
      nbImpCases <- rpois(1, rate.import.case)
      if(nbImpCases>0){
        ## dates of imported cases
        res$onset <- c(res$onset, rep(t, nbImpCases))
        
        ## ancestries of the imported cases
        res$ances <- c(res$ances, rep(NA, nbImpCases))
        
        ## id of the imported cases
        newId <- seq(N+1, by=1, length=nbImpCases)
        res$id <- c(res$id, newId)
        
        ## status of new cases
        res$status[newId] <- "I"
        
        ## group of the imported cases
        res$group <- c(res$group, choose.group(nbImpCases))
        
        ## dna sequences of the new infections
        newSeq <- t(sapply(1:nbImpCases, function(i) seq.dupli(EVE, diverg.import)))
        res$dna <- rbind(res$dna, newSeq)
      }
      
      
      ## set recovered status ##
      res$status[res$id[(t-res$onset) >= t.clear]] <- "R"
      
      ## update nb of infected, recovered, etc.
      res$dynam$nrec[t+1] <- sum(res$status=="R")
      res$dynam$ninf[t+1] <- sum(res$status=="I")
      res$dynam$nsus[t+1] <- sum(res$status=="S")
    } # end for
    
    
    ## SHAPE AND RETURN OUTPUT ##
    ## data need to be reordered so that res$id is 1:res$n
    res$n <- nrow(res$dna)
    res$ances <- match(res$ances, res$id)
    res$id <- 1:res$n
    res$xy <- res$inf.xy # don't keep entire distribution, not right order anymore anyway
    res$inf.xy <- NULL # simpler to just call coords 'xy'
    res$status <- NULL # we don't need this
    res$recover <- t.clear+res$onset
    
    findNmut <- function(i){
      if(!is.na(res$ances[i]) && res$ances[i]>0){
        out <- dist.dna(res$dna[c(res$id[i],res$ances[i]),], model="raw")*ncol(res$dna)
      } else {
        out <- NA
      }
      return(out)
    }
    
    ##res$nmut <- sapply(1:res$n, function(i) dist.dna(res$dna[c(res$id[i],res$ances[i]),], model="raw"))*ncol(res$dna)
    res$nmut <- sapply(1:res$n, function(i) findNmut(i))
    res$ngen <- rep(1, length(res$ances)) # number of generations
    res$call <- match.call()
    ## if(tree){
    ##     res$tree <- fastme.ols(dist.dna(res$dna, model="TN93"))
    ##     res$tree <- root(res$tree,"1")
    ## }
    
    
    class(res) <- "simOutbreak"
    return(res)
    
  } # end simOutbreak
  as.igraph.simOutbreak <- function(x, edge.col="black", col.edge.by="dist", vertex.col="orange",
                                    edge.col.pal=NULL, annot=c("dist","n.gen"), sep="/", ...){
    ## if(!require(igraph)) stop("package igraph is required for this operation")
    ## if(!require(ape)) stop("package ape is required for this operation")
    if(!inherits(x,"simOutbreak")) stop("x is not a tTree object")
    ## if(!require(adegenet)) stop("adegenet is required")
    if(!col.edge.by %in% c("dist","n.gen","prob")) stop("unknown col.edge.by specified")
    
    ## GET DAG ##
    from <- as.character(x$ances)
    to <- as.character(x$id)
    isNotNA <- !is.na(from) & !is.na(to)
    vnames <- unique(c(from,to))
    vnames <- vnames[!is.na(vnames)]
    dat <- data.frame(from,to,stringsAsFactors=FALSE)[isNotNA,,drop=FALSE]
    out <- graph.data.frame(dat, directed=TRUE, vertices=data.frame(names=vnames))
    
    ## from <- as.character(x$ances)
    ## to <- as.character(x$id)
    ## dat <- data.frame(from,to,stringsAsFactors=FALSE)[!is.na(x$ances),]
    ## out <- graph.data.frame(dat, directed=TRUE)
    
    ## SET VERTICE INFO ##
    ## labels
    V(out)$label <- V(out)$name
    
    ## dates
    names(x$onset) <- x$id
    V(out)$date <- x$onset[V(out)$name]
    
    ## ## groups
    ## names(x$group) <- x$id
    ## V(out)$group <- x$group[V(out)$name]
    
    ## colors
    V(out)$color <- vertex.col
    ## V(out)$color <- fac2col(factor(V(out)$group), col.pal=vertex.col.pal)
    
    
    ## SET EDGE INFO ##
    ## genetic distances to ancestor
    E(out)$dist <- x$nmut[!is.na(x$ances)]
    
    ## number of generations to ancestor
    E(out)$ngen <- x$ngen[!is.na(x$ances)]
    
    ## colors
    if(is.null(edge.col.pal)){
      edge.col.pal <- function(n){
        return(grey(seq(0.75,0,length=n)))
      }
    }
    if(col.edge.by=="dist") edge.col <- num2col(E(out)$dist, col.pal=edge.col.pal, x.min=0, x.max=1)
    
    ## labels
    n.annot <- sum(annot %in% c("dist","n.gen"))
    lab <- ""
    if(!is.null(annot) && n.annot>0){
      if("dist" %in% annot) lab <- E(out)$dist
      if("n.gen" %in% annot) lab <- paste(lab, x$ngen[!is.na(x$ances)], sep=sep)
    }
    lab <- sub(paste("^",sep,sep=""),"",lab)
    E(out)$label <- lab
    
    ## SET LAYOUT ##
    attr(out, "layout") <- layout.fruchterman.reingold(out, params=list(minx=V(out)$date, maxx=V(out)$date))
    
    return(out)
  } # end as.igraph.simOutbreak
  
  
  ##GENERATE OUTBREAK
  #repeat outbreak if doesn't require minimum cases
  n <- 0; i <- 1
  N <- 10 # max times of generating outbreak
  w <- EpiEstim::discr_si(k = seq(0, outbreak_duration), 
                          mu = genint_params[1], 
                          sigma = genint_params[2])
  while(n < minimum_case){
    myoutbreak <- tryCatch(theoutbreak(R0 = R0,
                                       infec.curve = w,
                                       n.hosts = initial_sus,
                                       duration = outbreak_duration,
                                       mu.transi = mutation_rate),
                           error = function(e){
                             list(n = -1, msg = conditionMessage(e))
                           })
    n <- myoutbreak$n
    if(n < minimum_case){
      if(n == -1) cat("Attempt", i, ":", "error detected.", myoutbreak$msg, "\n")
      else cat("Attempt", i, ":", "Outbreak doesn't exceed minimum cases.", "\n")
      cat("Repeating simulation.", "\n")
    } else cat("Attempt", i, ":", "Generate outbreak with ncases =", n, "\n")
    cat("\n")
    if(i == N){
      cat("Try to change smaller the minimum_case.")
      break
    }
    i <- i+1
  }
  
  ##GET THE OUTPUT: EPIDATA, WIWDATA, DNA
  epidata <- data.frame(inf.ID = paste("ID", myoutbreak$id, sep = "_"),
                        inf.times = myoutbreak$onset,
                        rec.times = myoutbreak$recover,
                        inf.source = paste("ID", ifelse(is.na(myoutbreak$ances), 0, myoutbreak$ances), sep = "_"))
  
  wiwdata <- epidata %>%
    select(inf.source) %>%
    left_join((epidata %>% select(inf.ID, inf.times)), by = c("inf.source" = "inf.ID")) %>%
    dplyr::rename(source.times = inf.times) %>%
    bind_cols((epidata %>% select(-c("inf.source")))) %>%
    mutate(si = inf.times-source.times) %>%
    select(inf.source, inf.ID, si)
  
  aligndata <- myoutbreak$dna
  row.names(aligndata) <- epidata$inf.ID
  aligndata <- as.list(aligndata)
  
  
  #add genomic distance to wiw data
  wiwdata$distance <- NA
  x <- ape::dist.dna(aligndata, as.matrix = TRUE)
  for(i in 1:nrow(wiwdata)){
    if(!is.na(wiwdata$si[i])){
      wiwdata$distance[i] <- x[wiwdata$inf.source[i], wiwdata$inf.ID[i]] 
    }
  }
  
  ##EXPORT DATAFRAME & DNABIN
  if(.export){
    write.csv(epidata, file = "sim_epidata.csv", row.names = F)
    write.csv(wiwdata, file = "sim_wiwdata.csv", row.names = F)
    ape::write.FASTA(aligndata, file = "sim_aligndata.fasta")
  }
  
  ##PLOT THE OUTBREAK
  if(.plot){
    gg <- as.igraph.simOutbreak(myoutbreak)
    par(mfrow = c(1,2), mar=c(0,0,0,0))
    plot(gg, layout=layout.reingold.tilford(gg), edge.arrow.size = .1,
         edge.arrow.width = .5, vertex.label = NA, edge.label.cex = .5,
         vertex.size = 2, edge.label = E(gg)$dist)
    plot(gg, layout=layout_with_fr(gg), edge.arrow.size = .1,
         edge.arrow.width = .5, vertex.label.cex = .3, edge.label.cex = .5,
         vertex.size = 3, vertex.label = NA, edge.label = wiwdata$si[!is.na(wiwdata$si)])
    par(mfrow = c(1,1))
  }
  
  ##GET THE OUTPUT
  output <- list(epidata = epidata, wiwdata = wiwdata, aligndata = aligndata)
  class(output) <- "outbreak"
  return(output)
}
asOutbreak <- function(case_id,
                       infection_time,
                       alignment = NA){
  aligndata <- alignment
  epidata <- data.frame(inf.ID = case_id,
                        inf.times = infection_time) %>%
    filter(inf.ID %in% names(aligndata))
  out <- list(epidata = epidata, aligndata = aligndata)
  class(out) <- "outbreak"
  return(out)
}
createTransCloud <- function(outbreak,
                             onset_interval_diff = c(0, Inf),
                             max_genetic_dist = 1e-3,
                             dnadist_model = "TN93"){
  if(class(outbreak)!="outbreak") stop("outbreak must be of class outbreak")
  
  myepidata <- outbreak$epidata
  myaligndata <- outbreak$aligndata
  pairs <- transtreesampler::id_transmission_pairs(data = myepidata,
                                                   aln = myaligndata,
                                                   max_genetic_dist = max_genetic_dist,
                                                   min_day_diff = onset_interval_diff[1],
                                                   onset_date_col = "inf.times",
                                                   case_id_col = "inf.ID",
                                                   max_day_diff = onset_interval_diff[2],
                                                   model = dnadist_model)
  pairs$onset_diff <- as.numeric(pairs$onset_diff)
  if(is.infinite(onset_interval_diff[2])) {
    maxinterval <- max(pairs$onset_diff) + 1
  } else maxinterval <- onset_interval_diff[2]
  if(is.infinite(max_genetic_dist)) {
    maxgendist <- max(pairs$distance) + .1
  } else maxgendist <- max_genetic_dist
  #maxgendist <- max_genetic_dist
  wtfn <- function(distance, onset_diff){
    if(maxgendist == 0 & maxinterval == 0) wt <- 1
    else if(maxgendist == 0) wt <- (maxinterval-onset_diff)/maxinterval
    else if(maxinterval == 0) wt <- (maxgendist-distance)/maxgendist
    else wt <- (maxinterval-onset_diff)/maxinterval + (maxgendist-distance)/maxgendist
    return(wt)
  }
  
  pairs %>%
    group_by(case_j) %>%
    mutate(weight = wtfn(distance, onset_diff),
           prob = weight/sum(weight)) %>%
    ungroup() %>%
    select(-c("weight"))
}

sampleTransTree <- function(transcloud,
                            count){
  ncore <- parallel::detectCores()
  cl <- makeCluster(ncore, type="SOCK")
  registerDoSNOW(cl)
  trees <- foreach(i = 1:count,
                   .packages = c("dplyr")) %dopar% {
                     transcloud %>%
                       group_by(case_j) %>%
                       slice_sample(n = 1, weight_by = prob) %>%
                       pull(onset_diff)
                   }
  stopCluster(cl)
  return(trees)
}

wiwAdj <- function(epidata, p){
  cases <- sources <- tcases <- M <- si <- c()
  epidata$inf.source[is.na(epidata$inf.source)] <- "ID_0"
  
  for(cs in epidata$inf.ID){
    u <- rbinom(1, 1, p)
    if(u == 0) next
    cases <- c(cases, cs)
    tcs <- with(epidata, inf.times[inf.ID == cs])
    tcases <- c(tcases, tcs)
    inf <- with(epidata, inf.source[inf.ID == cs])
    if(inf == "ID_0") u <- 1
    else u <- rbinom(1, 1, p)
    
    k <- 0
    while(u == 0){
      k <- k+1
      newcs <- inf
      inf <- with(epidata, inf.source[inf.ID == newcs])
      if(length(inf) == 0 || inf == "ID_0"){
        inf <- "ID_0"
        u <- 1
      } else u <- rbinom(1, 1, p)
    }
    tinf <- with(epidata, inf.times[inf.ID == inf])
    if(length(tinf) == 0) tinf <- NA
    sources <- c(sources, inf)
    si <- c(si, tcs-tinf)
    M <- c(M, k)
  }
  df <- data.frame(inf.ID = cases, inf.source = sources, inf.times = tcases, M = M, si = si)
  return(df)
}
simContacts <- function(nSample, par, R0, .plot = F){ 
  # generate observed serial interval with par
  x1 <- rcgg(round(par[4]*nSample), par[1], par[2], par[3])
  x2 <- rfgd(round((1-par[4])*nSample), par[1], par[2])
  si <- c(x1, x2)
  #si <- c(x1, x2)
  #
  
  # generate contacts
  case_i <- case_j <- c()
  pools <- bins <- c()
  ti <- tj <- c()
  nI <- 0 # acrually at t=0, we have ID_0, this is for indexing only
  nSus <- nSample
  #
  
  # generate and sample descendant of ID_0
  a <- max(round(rpois(1, R0) * par[3]), 1)
  case_i <- c(case_i, rep("ID_0", a))
  case_j <- c(case_j, paste("ID", (nI+1):(nI+a), sep = "_"))
  ti <- rep(0, a)
  tj <- si[1:a]
  pools <- case_j
  #
  
  # update dynamic
  nI <- nI + a
  nSus <- nSus - a
  #
  
  while(nSus > 0){
    if(length(pools) == 0){
      nI <- nI + 1
      inf <- paste("ID", nI, sep = "_") # import case
      a <- max(round(rpois(1, R0) * par[3]), 1)
    } else{
      inf <- sample(pools, 1)
      a <- round(rpois(1, R0) * par[3])
    } 
    
    # update case i and j
    pools <- pools[pools != inf]
    
    if(nSus < a) a <- nSus # to make sure we can only sample nSample cases
    
    #generate & sample decendent of inf
    if(a > 0){
      case_j <- c(case_j, paste("ID", (nI+1):(nI+a), sep = "_"))
      case_i <- c(case_i, rep(inf, a))
      
      if(inf %in% case_j) tmp <- tj[which(case_j == inf)]
      else tmp <- abs(max(tj)-abs(rnorm(1)))
      ti <- c(ti, rep(tmp, a))
      tj <- c(tj, tmp + si[which(case_i == inf)])
      
      pools <- c(pools, paste("ID", (nI+1):(nI+a), sep = "_"))
      nI <- nI + a
      nSus <- nSus - a
    }
  }
  tc <- data.frame(case_i = case_i, case_j = case_j, ti = ti, tj = tj, onset_diff = si, prob = 1)
  epi <- data.frame(inf.ID = case_j, inf.source = case_i, inf.times = tj)
  #
  
  # draw graph
  if(.plot){
    epi.nodes <- data.frame(id = unique(c(case_i, case_j))) %>% 
      left_join((tc %>% group_by(case_i) %>% summarise(value = n())), by = c("id"="case_i")) %>%
      mutate(value = ifelse(is.na(value), 1, value+1), font.size = 0, color = "coral")
    epi.edges <- data.frame(from = tc$case_i,
                            to = tc$case_j,
                            length = tc$onset_diff+5,
                            arrows = "to",
                            color = "black")
    graph <- list(nodes = epi.nodes, edges = epi.edges)
  } else graph <- NULL
  
  return(list(epi = epi, tt = tc, graph = graph))
}
#
#



