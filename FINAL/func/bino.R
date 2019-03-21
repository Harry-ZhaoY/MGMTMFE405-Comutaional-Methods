# combination
nCr <- function(n, r){
  if(r==0){
    return(1)
  } else {
    return(prod(sapply(1:r, function(i)(n-r+i)/(i))))
  }
}
# call option payoffs through binomial tree
binoCall <- function(u, d, p, s0, x, rf, t, n){
  dt <- t/n
  payoffs <- vector()
  for (k in 0:n){
    st <- s0*u^k*d^(n-k)
    payoffs[k+1] <- nCr(n,k)*p^k*(1-p)^(n-k)*ifelse(st>x, st-x, 0)
  }
  call <- sum(payoffs)*exp(-rf*t)
  return(call)
}
# put option binomial model
putBino <- function(u, d, p, s0, x, rf, t, n, American){
  dt <- t/n
  if (American){
    put <- rep(0,(n+2))
    for (step in 0:n){
      k <- (n-step):0
      st <- s0*u^k*d^((n-step)-k)
      cv <- exp(-rf*dt)*na.omit(c(NA,put)*p + c(put,NA)*(1-p))
      put <- ifelse(cv>x-st, cv, ifelse(st<x, x-st, 0))
    }
  } else {
    payoffs <- vector()
    for (k in 0:n){
      st <- s0*u^k*d^(n-k)
      payoffs[k+1] <- nCr(n,k)*p^k*(1-p)^(n-k)*ifelse(st<x, x-st, 0)
    }
    put <- sum(payoffs)*exp(-rf*t)
    
  }
  return(put)
}
# trinomial call
trinoCall <- function(u, d, pU, pD, s0, x, rf, t, n, logPrice){
  powU <- c(n:0, rep(0,n))
  powD <- c(rep(0,(n)), 0:(n))
  if(logPrice){
    st <- exp(log(s0)+u*powU+d*powD)
  } else {
    st <- s0*u^powU*d^powD
  }
  call <- ifelse(st>x, st-x, 0)
  for (step in 1:n){
    call <- exp(-rf*dt)*na.omit((c(NA, NA, call)*pU+c(NA, call, NA)*(1-pU-pD)+c(call, NA, NA)*pD))
  }
  return(call)
}

