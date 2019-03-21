defaultOption <- function(paths, lambda2, t, L0, r0, delta, alpha, epsilon){
  n <- ncol(paths)
  m <- nrow(paths)
  mons <- t*12
  # find r
  r <- (r0+delta*lambda2)/12
  # find pmt
  pmt <- (L0*r)/(1-(1/(1+r))^mons)
  # find Lt
  a <- pmt/r
  b <- pmt/(r*(1+r)^mons)
  c <- 1+r
  Lt_mon <- a-b*c^(0:mons)
  # find qt
  beta <- (epsilon-alpha)/t
  qt_mon <- alpha + beta*(0:mons)
  # use linear interpolation for continuous loan values
  Lt <- approx(seq(0,n,n/mons), Lt_mon, n=n)[[2]]
  qt <- approx(seq(0,n,n/mons), qt_mon, n=n)[[2]]
  # find Nt
  Nt <- matrix(rpois(n*m, lambda2), m, n, byrow = T)
  # find Q
  Q <- unlist(apply(paths, 1, function(x) which(x<=qt*Lt)[1]))
  # find S
  S <- unlist(apply(Nt, 1, function(x) which(x>0)[1]))
  # find value of default option, default prob, and avg default time
  option <- vector()
  tao <- vector()
  for (i in 1:m){
    if (!is.na(Q[i])&!is.na(S[i])){
      # if Q happens before S
      if (Q[i]<=S[i]){
        option[i] <- max(Lt[Q[i]]-epsilon*paths[i,Q[i]], 0)
        tao[i] <- Q[i]
      }
      # if S happens before Q
      else{
        option[i] <- abs(Lt[S[i]]-epsilon*paths[i,S[i]])
        tao[i] <- S[i]
      }
    }
    # if only Q happens
    else if (!is.na(Q[i])){
      option[i] <- max(Lt[Q[i]]-epsilon*paths[i,Q[i]], 0)
      tao[i] <- Q[i]
    }
    # if only S happens
    else if (!is.na(S[i])){
      option[i] <- abs(Lt[S[i]]-epsilon*paths[i,S[i]])
      tao[i] <- S[i]
    }
    # if no default happens
    else {
      option[i] <- 0
      tao[i] <- NA
    }
  }
  prob <- sum(option!=0)/m
  # return
  return(c(value = mean(option), prob = prob, tao =mean(tao, na.rm = T)))
}
