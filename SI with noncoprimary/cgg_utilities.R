dcgg <- function(x,a,b,pi){
  #setup
  if(a<0 || b<0 || pi<0){
    par <- c(a,b,pi)
    par[which(par<0)] <- 0
  }
  if(pi>1) pi <- 1
  #Find maximum m within tol=1e-10
  k <- 0; tol <- 1e-10
  y <- 1/gamma(a)
  v <- b^a * x^a * (1-pi)
  while(y>tol && y<Inf){
    k <- k+1
    y <- tryCatch(v^k/gamma((k+1)*a), error=function(e) return(0), warning=function(w) return(0))
    if(is.nan(y)) y <- 0 #break the loop
    #NaN happens when gamma((m+1)*a)=Inf
  }
  if(k>1) mmax <- k-1
  else mmax <- 0
  
  f <- 0
  for(kk in 0:mmax) f <- f + dgamma(x,(kk+1)*a,b)*dgeom(kk,pi)
  f[f==0] <- 1e-300 #to avoid log(0)
  return(f)
}

#log-likelihood 
llcgg <- function(params, data){
  
  # switch to shape alpha, rate beta
  alpha <- (params[1]^2)/(params[2]^2)
  beta <- params[1]/(params[2]^2)
  params[c(1,2)] <- c(alpha, beta)
  
  #likelihood
  f <- numeric(length(data))
  for(i in 1:length(data)) f[i] <- dcgg(data[i], params[1], params[2], params[3])
  logf <- sum(log(f))
  
  return(logf)
}

# log-posterior
lpcgg <- function(params, data, pi_prior=NULL){
  
  if (is.null(pi_prior)){stop("For the chosen method, you must provide a prior for pi, named 'pi_prior'")}
  
  # switch to shape alpha, rate beta
  alpha <- (params[1]^2)/(params[2]^2)
  beta <- params[1]/(params[2]^2)
  params[c(1,2)] <- c(alpha, beta)
  
  #likelihood
  f <- numeric(length(data))
  for(i in 1:length(data)) f[i] <- dcgg(data[i], params[1], params[2], params[3])
  logf <- sum(log(f))
  #prior
  prior <- (dbeta(params[3],pi_prior[1],pi_prior[2]))
  lprior <- ifelse(prior==0, log(1e-300), log(prior))
  return(logf+lprior)
}


#call utilities.R
siCGG <- function(si, init=list(mu,sigma,pi=NULL), lower=c(1e-5,1e-5,1e-5), 
                  upper=c(30,30,1), ci=NULL, pi_prior=NULL){
    #setup
    if(is.null(ci)) ci <- 0.95
    params0 <- c(init$mu, init$sigma)
    if(is.null(init$pi)) stop("Please provide initial value of pi.")
    params0 <- c(params0, init$pi)
    #end setup
    
    #estimate the parameters
    if (!is.null(pi_prior)){
      # prior on pi
      myestim <- tryCatch(optim(params0, lpcgg, data=si, pi_prior=pi_prior, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=T),
                          error=function(x){
                            e <- optim(params0, lpcgg, data=si, pi_prior=pi_prior, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=F)
                            return(e)
                          })
    }else{
        # pi is estimated
        myestim <- tryCatch(optim(params0, llcgg, data=si, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=T),
        error=function(x){
             e <- optim(params0, llcgg, data=si, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=F)
             return(e)
        })
    }
      
    ## params
    estim <- myestim$par
    names(estim) <- c("mu", "sigma", "pi")
    if(all(estim>lower) && all(estim<upper)) msg <- myestim$message
    else msg <- "ESTIMATES ARE EXACTLY AT THE BOUND."
  
    ##log-likelihood value
    logLik <- myestim$value
    ##hessian matrix
    hess <- myestim$hessian
    if(is.null(hess) || prod(diag(hess))==0){ 
      if (!is.null(pi_prior) || !is.null(w_prior)){ hess <- numDeriv::hessian(func=lpcgg, x=estim, data=si)
      } else{ hess <- numDeriv::hessian(func=llcgg, x=estim, data=si)}}
    varcov <- solve(-hess) #fisher information
    
    ##confidence interval
    z <- qnorm(ci+(1-ci)/2)
    matCI <- matrix(0, nrow=3, ncol=2)
    rownames(matCI) <- c("mu", "sigma", "pi")
    colnames(matCI) <- c("lower","upper")
    matCI[,1] <- estim-z*sqrt(abs(diag(varcov)))
    matCI[,2] <- estim+z*sqrt(abs(diag(varcov)))
    for(i in 1:nrow(matCI)){
        if(matCI[i,1]<0) matCI[i,1] <- 0
    }
    if(matCI[3,2]>1) matCI[3,2] <- 1
    
    ##output
    out <- list(estimates=estim, logLik=logLik, varcov=varcov, ci.matrix=matCI, message=msg)
    return(out)
}
