lsmc <- function(paths, strike, r, t, nSim, n, func, k){
  dt <- t/n
  # calculate exercise value at each time
  ev <- ifelse(strike-paths>0, strike-paths,0)
  # construct index
  ind <- cbind(diag(0, nSim, n-1), ifelse(ev[,n]>0, 1,0))
  # for loop through time
  for (i in 1:(n-2)){
    # define x
    x <- paths[,(n-i)]
    # define discount factors
    dsc <- exp(-r*dt*matrix(rep(1:i,each = nSim), nSim, i))
    # define y
    y <- apply(matrix(ind[,(n-i+1):n]*ev[,(n-i+1):n]*dsc,
                      nSim, i),1,sum)
    # define basis functions
    if (func == "Laguerre"){
      L <- function(x, k){
        out <- switch(k,
                      exp(-x/2),
                      exp(-x/2)*(1-x),
                      exp(-x/2)*(1-2*x+x^2/2),
                      exp(-x/2)*(1-3*x+(3*x^2)/2-x^3/6))
        return(out)
      }
    } 
    else if (func == "Hermite"){
      L <- function(x, k){
        out <- switch(k,
                      x^0,
                      2*x,
                      4*x^2-2,
                      8*x^3-12*x)
        return(out)
      }
    } 
    else if (func == "Monomials"){
      L <- function(x, k){
        return(x^(k-1))
      }
    } 
    else {
      stop("Incorrect function type!")
    }
    # Find matrix A and vector b
    A <- diag(0, k, k)
    b <- vector()
    for (p in 1:k){
      for (q in 1:k){
        A[p,q] <- sum(L(x,p)*L(x,q))
      }
      b[p] <- sum(y*L(x,p))
    }
    # find a_i coefficients, using cholesky roots to inverse matrix
    a <- chol2inv(chol(A))%*%b
    
    # find expected continuation values
    est <- vector()
    for (p in 1:k){
      est <- cbind(est, L(x,p))
    }
    expCv <- est%*%a
    
    # rewirte index
    ind[,(n-i)] <- ifelse(ev[,(n-i)]>expCv, 1, 0)
    # write off 1s in previous
    ind <- t(apply(ind, 1, function(x) {
      if (sum(x)>1){
        replace(x, (which(x>0)[1]+1):n, 0)
      } 
      else {x}
    }))
  } 
  # create discount matrix
  dsct <- exp(-r*dt*matrix(rep(0:(n-1),each = nSim), nSim, n))
  
  # compute option value at time 0
  v0 <- mean(apply(ind*dsct*ev, 1, sum))
  return(v0)
}