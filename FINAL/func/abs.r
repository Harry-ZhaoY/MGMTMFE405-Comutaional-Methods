zcBond_cir <- function(r0, sig, kappa, r_bar, t){
  h1 <- sqrt(kappa^2+2*sig^2)
  h2 <- (kappa+h1)/2
  h3 <- (2*kappa*r_bar)/sig^2
  A <- ((h1*exp(h2*t))/(h2*(exp(h1*t)-1)+h1))^h3
  B <- (exp(h1*t)-1)/(h2*(exp(h1*t)-1)+h1)
  return(A*exp(-B*r0))
}

cirPath <- function(r0, sig, kappa, r_bar, t, nSim){
  dt <- 1/12 
  n <- t/dt+1
  r_sim <- matrix(rep(r0, nSim), ncol = 1)
  for (s in 2:n){
    dw <- rnorm(nSim)*sqrt(dt)
    # use partial truncation
    newCol <- r_sim[,s-1] + kappa*(r_bar - r_sim[,s-1])*dt + sig*sqrt(ifelse(r_sim[,s-1]>0,r_sim[,s-1],0))*dw
    r_sim <- cbind(r_sim, newCol)
  }
  return(r_sim)  
}

findMBS <- function(pv0=100000, r0=.078, r_bar=.08, wac =.08, kappa=.6, sig=.12, years=30, nSim=30000){
  # function to find discount factor
  disFac <- function(r, t, n){
    factor <- 1/(1-(1+r)^(-n+(t-1)))
    return(factor)
  }
  
  # set constants
  r <- wac/12
  N <- years*12+1
  # find CPR_t related terms
  t_seq <- seq(0, years*12)
  SG_t <- ifelse(t_seq/30<1, t_seq/30, 1)
  SY_t <- rep(c(.94,.76,.74,.95,.98,.92,.98,1.1,1.18,1.22,1.23,.98), years+1)
  
  # generate interest rate paths
  r_t <- cirPath(r0, sig, kappa, r_bar, years, nSim)
  # set containers for PV_t and PV_CF_t
  PV_t <- cbind(rep(pv0,nSim), diag(0, nSim, N))
  PV_CF_t <- diag(0, nSim, N)
  for (s in 1:N){
    # find r_t-1(10)
    r10 <- -log(zcBond_cir(r_t[,s], sig, kappa, r_bar, 10))/10
    RI_s <- 0.28+0.14*atan(-8.57 + 430*(wac - r10))
    BU_s <- 0.3+0.7*(PV_t[,s]/pv0)
    # find CPR_s
    CPR_s <- RI_s*BU_s*SG_t[s]*SY_t[s]
    # find cahs flow s
    SP_s <- PV_t[,s]*r*(disFac(r, s, N)-1)
    IP_s <- PV_t[,s]*r
    TPP_s <- SP_s + (PV_t[,s] - SP_s) * (1-(1-CPR_s)^(1/12))
    # update PV_t
    PV_t[,s+1] <- PV_t[,s] - TPP_s
    # update PV_CF_t
    r_0_s <- if (s == 1) matrix(r_t[,1], ncol = 1) else  r_t[,1:s]
    PV_CF_t[,s] <- (TPP_s + IP_s)*exp(-apply(r_0_s*(1/12), 1, sum))
  }
  P_0 <- mean(apply(PV_CF_t, 1, sum))
  return(P_0)
}

fitOAS <- function(x, r_t, pv0=100000, wac =.08, expPv = 110000){
  # function to find discount factor
  disFac <- function(r, t, n){
    factor <- 1/(1-(1+r)^(-n+(t-1)))
    return(factor)
  }
  
  # set constants
  r <- wac/12
  N <- (ncol(r_t)-1)+1
  nSim <- nrow(r_t)
  
  # find CPR_t related terms
  t_seq <- seq(0, N-1)
  SG_t <- ifelse(t_seq/30<1, t_seq/30, 1)
  SY_t <- rep(c(.94,.76,.74,.95,.98,.92,.98,1.1,1.18,1.22,1.23,.98), (N-1)/12+1)
  
  # add spread
  r_t <- r_t+x
  # set containers for PV_t and PV_CF_t
  PV_t <- cbind(rep(pv0,nSim), diag(0, nSim, N))
  PV_CF_t <- diag(0, nSim, N)
  for (s in 1:N){
    # find r_t-1(10)
    r10 <- -log(zcBond_cir(r_t[,s], 0.12, 0.6, 0.08, 10))/10
    RI_s <- 0.28+0.14*atan(-8.57 + 430*(wac - r10))
    BU_s <- 0.3+0.7*(PV_t[,s]/pv0)
    # find CPR_s
    CPR_s <- RI_s*BU_s*SG_t[s]*SY_t[s]
    # find cahs flow s
    SP_s <- PV_t[,s]*r*(disFac(r, s, N)-1)
    IP_s <- PV_t[,s]*r
    TPP_s <- SP_s + (PV_t[,s] - SP_s) * (1-(1-CPR_s)^(1/12))
    # update PV_t
    PV_t[,s+1] <- PV_t[,s] - TPP_s
    # update PV_CF_t
    r_0_s <- if (s == 1) matrix(r_t[,1], ncol = 1) else  r_t[,1:s]
    PV_CF_t[,s] <- (TPP_s + IP_s)*exp(-apply(r_0_s*(1/12), 1, sum))
  }
  
  P_0 <- mean(apply(PV_CF_t, 1, sum))
  return(P_0-expPv)
}