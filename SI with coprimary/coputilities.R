reltol <- function(past, now) abs(past-now)/past
rcgg <- function(n,a,b,pi){
  if(a<0 || b<=0 || pi<=0 || pi>1){
    if(a<0) stop("Invalid parameter: a<0.")
    else if(b<=0) stop("Invalid parameter: b<=0")
    else stop("Invalid parameter: p not in (0,1].")
  }
  m <- rgeom(n,pi)
  x <- rgamma(n, (m+1)*a, b)
  return(x)
}
rgdd <- function(n,a,b){
  if(length(a)!=2 || length(b)!=2) stop("Provide a=(a1,a2) and b=(b1,b2).")
  x1 <- rgamma(n, a[1], b[1])
  x2 <- rgamma(n, a[2], b[2])
  return(x1-x2)
}
rfgd <- function(n,a,b){
  if(a<=0) stop("May need to consider choosing a>0.")
  if(b<=0) stop("May need to consider choosing b>0.")
  x1 <- rgamma(n, a, b)
  x2 <- rgamma(n, a, b)
  return(abs(x1-x2))
}
randomSI <- function(n,a,b,pi,w){
  x <- numeric(n)
  z <- rbinom(n, 1, w)
  n1 <- which(z==1)
  n2 <- which(z==0)
  x[n1] <- rcgg(length(n1), a, b, pi)
  x[n2] <- rfgd(length(n2), a, b)
  return(list(si=x, z=z))
}
dcgg2 <- function(x,a,b,pi){
  #setup
  if(a<=0) stop("NaN produced for a<=0.")
  if(b<=0) stop("Provide b>0.")
  if(pi<=0 || pi>1) stop("Provide pi in (0,1].")
  
  gt <- function(xi,a,b,pi){
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
    if(length(which(is.na(f)))>1){
      cat("There's NaN when computing dcgg2 at x =", xi, "\n")
      cat("print(par):", a, b, pi, "\n")
      stop("")
    }
    return(sum(f))
  }
  ft <- tryCatch(sapply(1:length(x), function(i) gt(xi=x[i], a, b, pi)),
                 error=function(e){e})
  return(ft)
}
dgdd <- function(x,a,b){
  if(length(a)!=2 || length(b)!=2) stop("Provide a=(a1,a2) and b=(b1,b2).")
  if(a[1]<=0 || a[2]<=0) stop("NaN produced for a<=0.")
  if(b[1]<=0 || b[2]<=0) stop("Provide b>0.")
  gt <- function(ti,a,b){
    cons <- b[1]^a[1] * b[2]^a[2] / (gamma(a[1])*gamma(a[2]))
    if(ti>=0){
      integrand <- function(x,s) x^(a[1]-1) * (x-s)^(a[2]-1) * exp(-sum(b)*x)
      int <- tryCatch(integrate(integrand, lower=ti, upper=Inf, s=ti)$value, error=function(e) return(0))
          #when the integral is divergent, set the density=0
      out <- cons*exp(b[2]*ti)*int
    } 
    else{
      integrand <- function(x,s) x^(a[2]-1) * (x+s)^(a[1]-1) * exp(-sum(b)*x)
      int <- tryCatch(integrate(integrand, lower=-ti, upper=Inf, s=ti)$value, error=function(e) return(0))
          #when the integral is divergent, set the density=0
      out <- cons*exp(-b[1]*ti)*int
    } 
    if(is.nan(out)) out <- 0 #happens because we may multiply 0 and Inf
    return(out)
  }
  ft <- tryCatch(sapply(1:length(x), function(i) gt(ti=x[i], a, b)),
                 error=function(e){
                   cat("Problem when evaluation x at:", x[i], "\n")
                   cat("print(par):", a, b, "\n")
                   msg <- conditionMessage(e)
                   return(msg)
                   })
  if(class(ft)!="numeric") stop("")
  return(ft)
}
dfgd <- function(x,a,b){
  if(a<=0) stop("NaN produced for a<=0.")
  if(b<=0) stop("Provide b>0.")
  ft <- numeric(length(x))
  ft[which(x>=0)] <- dgdd(x[x>=0], c(a,a), c(b,b))
  return(2*ft)
}
dcop <- function(x,par,w){
  fx <- w*dcgg2(x, par[1], par[2], par[3]) + (1-w)*dfgd(x, par[1], par[2])
  return(fx)
}
pcgg2 <- function(x,a,b,pi){
  pt <- function(xi,a,b,pi){
    int <- integrate(f=dcgg2, lower=0, upper=xi, a=a, b=b, pi=pi)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], a, b, pi))
  Ft[which(Ft>1)] <- 1
  return(Ft)
}
pgdd <- function(x,a,b){
  pt <- function(xi,a,b){
    int <- integrate(f=dgdd, lower=-Inf, upper=xi, a=a, b=b)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], a, b))
  Ft[which(Ft>1)] <- 1
  return(Ft)
}
pfgd <- function(x,a,b){
  pt <- function(xi,a,b){
    int <- integrate(f=dfgd, lower=0, upper=xi, a=a, b=b)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], a, b))
  Ft[which(Ft>1)] <- 1
  return(Ft)
}
pcop <- function(x,par,w){
  pt <- function(xi,par,w){
    int <- integrate(f=dcop, lower=0, upper=xi, par=par, w=w)
    val <- int$value
    return(val)
  }
  Ft <- sapply(1:length(x), function(i) pt(x[i], a, b))
  Ft[Ft>1] <- 1
  return(Ft)
}

llOPTIM <- function(params, data){
  
  # Convert mu, sigma to alpha, beta
  alpha <- (params[1]^2)/(params[2]^2)
  beta <- params[1]/(params[2]^2)
  params[c(1,2)] <- c(alpha, beta)
  
  #setup
  if(params[1]<=0) params[1] <- 1e-300 #get error when computing gamma(0)
  if(params[2]<=0) params[2] <- 1e-300
  if(params[3]<=0) params[3] <- 1e-300
  if(params[3]>1) params[3] <- 1
  if(params[4]<0) params[4] <- 0
  if(params[4]>1) params[4] <- 1
      #this is only for optimization
  
  f1 <- dcgg2(data, params[1], params[2], params[3])
  f2 <- dfgd(data, params[1], params[2])
  f1[f1==0] <- 1e-300 #avoid log(0)
  f2[f2==0] <- 1e-300 #avoid log(0)
      #L-BFGS-B only allows finite obj. func.
  logf <- sum(log(params[4]*f1 + (1-params[4])*f2))
  #logf <- sum(params[4]*log(f1) + (1-params[4])*log(f2))
  return(logf)
}

postOPTIM <- function(params, data, pi_prior=NULL, w_prior=NULL){
  
  if (is.null(pi_prior) && is.null(w_prior)){warning('You chose method `prior`, but did not provide any priors')}
  
  # Convert mu, sigma to alpha, beta
  alpha <- (params[1]^2)/(params[2]^2)
  beta <- params[1]/(params[2]^2)
  params[c(1,2)] <- c(alpha, beta)
  
  #setup
  if(params[1]<=0) params[1] <- 1e-300 #get error when computing gamma(0)
  if(params[2]<=0) params[2] <- 1e-300
  if(params[3]<=0) params[3] <- 1e-300
  if(params[3]>1) params[3] <- 1
  if(params[4]<0) params[4] <- 0
  if(params[4]>1) params[4] <- 1
  #this is only for optimization
  
  f1 <- dcgg2(data, params[1], params[2], params[3])
  f2 <- dfgd(data, params[1], params[2])
  f1[f1==0] <- 1e-300 #avoid log(0)
  f2[f2==0] <- 1e-300 #avoid log(0)
  #L-BFGS-B only allows finite obj. func.
  logf <- sum(log(params[4]*f1 + (1-params[4])*f2))
  #logf <- sum(params[4]*log(f1) + (1-params[4])*log(f2))
  
  #prior
  if (!is.null(pi_prior)){
  prior <- (dbeta(params[3],pi_prior[1],pi_prior[2]))
  lprior <- ifelse(prior==0, log(1e-300), log(prior))
  }else{lprior <- 0}
  if (!is.null(w_prior)){
    priorw <- dbeta(params[4],w_prior[1], w_prior[2])
    lprior <- lprior + ifelse(priorw==0, log(1e-300), log(priorw))
  }

  return(logf+lprior)
}


OPTIMestim <- function(si,par0,lower,upper,ci, pi_prior=NULL, w_prior=NULL){
  #setup
  n <- length(si)
  param0 <- par0
  
  if (!is.null(pi_prior) || !is.null(w_prior)){
    myestim <- tryCatch(optim(par=param0, fn=postOPTIM, data=si, pi_prior=pi_prior, w_prior=w_prior, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=T), 
                      error=function(x){
                        e <- optim(par=param0, fn=postOPTIM, data=si, pi_prior=pi_prior, w_prior=w_prior, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=F)
                        return(e)
                      })
  }else{
    myestim <- tryCatch(optim(par=param0, fn=llOPTIM, data=si, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=T), 
                        error=function(x){
                          e <- optim(par=param0, fn=llOPTIM, data=si, method="L-BFGS-B", lower=lower, upper=upper, control=list("fnscale"=-1), hessian=F)
                          return(e)
                        })
    
  }
  
  param <- myestim$par
  names(param) <- c("mu", "sigma", "pi", "w")
  if(all(param>lower) && all(param<upper)) msg <- myestim$message
  else msg <- "ESTIMATES ARE EXACTLY AT THE BOUND."
  logLik <- myestim$value
  hess <- myestim$hessian
  if(is.null(hess) || prod(diag(hess))==0){ 
    if (!is.null(pi_prior) || !is.null(w_prior)){ hess <- numDeriv::hessian(func=postOPTIM, x=param, data=si)
                                          } else{hess <- numDeriv::hessian(func=llOPTIM, x=param, data=si)}}
  varcov <- solve(-hess) #fisher information
  se <- sqrt(abs(diag(varcov)))
  
  ##confidence interval
  z <- qnorm(ci+(1-ci)/2)
  matCI <- matrix(0, nrow=4, ncol=2)
  rownames(matCI) <- c("mu", "sigma", "pi", "w")
  colnames(matCI) <- c("lower","upper")
  matCI[,1] <- param-z*se
  matCI[,2] <- param+z*se
  matCI[,1][matCI[,1]<0] <- 0
  matCI[3:4,2][matCI[3:4,2]>1] <- 1
  
  out <- list(par=param, se=se, varcov=varcov, ci.matrix=matCI, logLik=logLik, message=msg)
  return(out)
}