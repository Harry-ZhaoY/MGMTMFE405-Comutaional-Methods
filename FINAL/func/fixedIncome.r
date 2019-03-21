# Zhao_Yanxiang_Projcet8
########################### FUNCTIONS ########################### 
##### BOND PRICING #####
# Zero-coupon bond
zcBond <- function(paths, fv){
  dt <- 1/360
  if (!is.null(ncol(paths))){
    pv <- fv*mean(exp(-apply(paths[,2:ncol(paths)]*dt, 1, sum)))
  } else{
    # deal with a single paths
    pv <- fv*mean(exp(sum(paths[2:length(paths)]*dt)))
  }
  return(pv)
}
# Coupon bond
cpnBond <- function(paths, pmtT, c, fv){
  pmtT <- pmtT*360
  pmtC <- c(rep(c, length(pmtT)-1), c+fv)
  pmtPV <- c()
  for (i in 1:length(pmtT)){
    pmtPV[i] <- zcBond(paths[,1:pmtT[i]],pmtC[i])
  }
  return(sum(pmtPV))
}
zcBond_cir <- function(r0, sig, kappa, r_bar, t, fv){
  h1 <- sqrt(kappa^2+2*sig^2)
  h2 <- (kappa+h1)/2
  h3 <- (2*kappa*r_bar)/sig^2
  A <- ((h1*exp(h2*t))/(h2*(exp(h1*t)-1)+h1))^h3
  B <- (exp(h1*t)-1)/(h2*(exp(h1*t)-1)+h1)
  return(fv*A*exp(-B*r0))
}
##### OPTION PRICING #####
# European option on ZC bond - vasicek
zcOption_vasicek <- function(paths, sig, kappa, r_bar, exerT, t, k, fv, type){
  dt <- 1/360
  exerT_day <- exerT*360+1
  rt <- paths[,exerT_day]
  B <- 1/kappa*(1-exp(-kappa*(t-exerT)))
  A <- exp((r_bar - sig^2/(2*kappa^2))*(B-(t-exerT)) - sig^2/(4*kappa)*B^2)
  bond_exerT <- fv*(A*exp(-B*rt))
  payoff <- switch (type,
    "call" = ifelse(bond_exerT>k, bond_exerT-k, 0),
    "put"  = ifelse(bond_exerT<k, k-bond_exerT, 0)
  )
  option <- mean(payoff*exp(-apply(paths[,2:exerT_day]*dt, 1, sum)))
  return(option)
}
# European call option on ZC bond - CIR explicit
zcCall_cir <- function(r0, sig, kappa, r_bar, exerT, t, k, fv){
  h1 <- sqrt(kappa^2+2*sig^2)
  h2 <- (kappa+h1)/2
  h3 <- (2*kappa*r_bar)/sig^2
  A_TS <- ((h1*exp(h2*(t-exerT)))/(h2*(exp(h1*(t-exerT))-1)+h1))^h3
  B_TS <- (exp(h1*(t-exerT))-1)/(h2*(exp(h1*(t-exerT))-1)+h1)
  theta <- sqrt(kappa^2+2*sig^2)
  phi <- 2*theta/(sig^2*(exp(theta*t)-1))
  psi <- (kappa + theta)/sig^2
  r_star <- log(A_TS/(k/fv))/B_TS
  P_tS <- zcBond_cir(r0, sig, kappa, r_bar, t, fv)/fv
  P_tT <- zcBond_cir(r0, sig, kappa, r_bar, exerT, fv)/fv
  # find the call option price with chi-dist
  x1 <- 2*r_star*(phi + psi + B_TS)
  q1 <- (2*phi^2*r0*exp(theta*t))/(phi+psi+B_TS)
  p <- (4*kappa*r_bar)/sig^2
  x2 <- 2*r_star*(phi + psi)
  q2 <- (2*phi^2*r0*exp(theta*t))/(phi+psi)
  # call option price
  C <- fv*P_tS*pchisq(x1, p, q1)-k*P_tT*pchisq(x2, p, q2)
  return(C)
}
# European option on ZC bond
zcOption <- function(paths, sig, kappa, r_bar, exerT, t, k, fv, type, method, mSim){
  rt <- paths[,ncol(paths)]
  paths_before <- paths
  option <- c()
  for (i in 1:length(rt)){
    paths_after <- switch (method,
                           "Vasicek" = vasicekPath(rt[i], sig, kappa, r_bar, t-exerT, mSim),
                           "CIR" = cirPath(rt[i], sig, kappa, r_bar, t-exerT, mSim)
    )
    bond_exerT <- zcBond(paths_after, fv)
    payoff <- switch (type,
                      "call" = ifelse(bond_exerT>k, bond_exerT-k, 0),
                      "put"  = ifelse(bond_exerT<k, k-bond_exerT, 0)
    )
    option[i] <- zcBond(paths_before[i,], payoff)
  }
  return(mean(option))
}
# European option on coupon bond
cpnOption <- function(paths, sig, kappa, r_bar, exerT, t, k, c, fv, type, method, mSim){
  rt <- paths[,ncol(paths)]
  paths_before <- paths
  option <- c()
  for (i in 1:length(rt)){
    paths_after <- switch (method,
      "Vasicek" = vasicekPath(rt[i], sig, kappa, r_bar, t-exerT, mSim),
      "CIR" = cirPath(rt[i], sig, kappa, r_bar, t-exerT, mSim)
    )
    pmtT <- seq(0.5, t, 0.5)-exerT # so far only for exerT < 0.5
    bond_exerT <- cpnBond(paths_after, pmtT, c, fv)
    payoff <- switch (type,
      "call" = ifelse(bond_exerT>k, bond_exerT-k, 0),
      "put"  = ifelse(bond_exerT<k, k-bond_exerT, 0)
    )
    option[i] <- zcBond(paths_before[i,], payoff)
  }
  return(mean(option))
}
# European opion on zc bond - G2++
zcOption_gcpp <- function(paths, x0, y0, r0, phi_t, rho, a, b, sig, eta, exerT, t, k, fv, type, mSim){
  rt <- paths[,ncol(paths)]
  paths_before <- paths
  option <- c()
  for (i in 1:length(rt)){
    paths_after <- g2ppPath(x0, y0, rt[i], phi_t, rho, a, b, sig, eta, t, mSim)
    bond_exerT <- zcBond(paths_after, fv)
    payoff <- switch (type,
                      "call" = ifelse(bond_exerT>k, bond_exerT-k, 0),
                      "put"  = ifelse(bond_exerT<k, k-bond_exerT, 0)
    )
    option[i] <- zcBond(paths_before[i,], payoff)
  }
  return(mean(option))
}