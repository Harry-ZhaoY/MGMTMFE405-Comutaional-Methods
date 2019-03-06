# Zhao_Yanxiang_Projcet8
########################### FUNCTIONS ########################### 
##### PATH #####
# Vasicek model
vasicekPath <- function(r0, sig, kappa, r_bar, t, nSim){
  dt <- 1/360 # assume 30/360 convension
  n <- t/dt
  sim <- list()
  for (i in 1:(nSim/2)){
    dw <- sqrt(dt)*rnorm(n)
    r_sim1 <- rep(r0,n)
    r_sim2 <- rep(r0,n)
    for (s in 2:n){
      # create anthithetic paths
      r_sim1[s] <- r_sim1[s-1] + kappa*(r_bar - r_sim1[s-1])*dt + sig*dw[s]
      r_sim2[s] <- r_sim2[s-1] + kappa*(r_bar - r_sim2[s-1])*dt - sig*dw[s]
    }
    sim_i <- list(r_sim1, r_sim2)
    sim <- append(sim, sim_i)
  }
  out <- matrix(unlist(sim), nSim, n, byrow = T)
  return(out)  
}
# CIR model
cirPath <- function(r0, sig, kappa, r_bar, t, nSim){
  dt <- 1/360 # assume 30/360 convension
  n <- t/dt
  sim <- list()
  for (i in 1:(nSim/2)){
    dw <- sqrt(dt)*rnorm(n)
    r_sim1 <- rep(r0,n)
    r_sim2 <- rep(r0,n)
    for (s in 2:n){
      # create anthithetic paths
      # use partial truncation
      r_sim1[s] <- r_sim1[s-1] + kappa*(r_bar - r_sim1[s-1])*dt + sig*sqrt(ifelse(r_sim1[s-1]>0,r_sim1[s-1],0))*dw[s]
      r_sim2[s] <- r_sim1[s-1] + kappa*(r_bar - r_sim2[s-1])*dt - sig*sqrt(ifelse(r_sim2[s-1]>0,r_sim2[s-1],0))*dw[s]
    }
    sim_i <- list(r_sim1, r_sim2)
    sim <- append(sim, sim_i)
  }
  out <- matrix(unlist(sim), nSim, n, byrow = T)
  return(out)  
}
# G2++ model
g2ppPath <- function(x0, y0, r0, phi_t, rho, a, b, sig, eta, t, nSim){
  dt <- 1/360 # assume 30/360 convension
  n <- t/dt
  sim <- list()
  for (i in 1:nSim){
    dw1 <- sqrt(dt)*rnorm(n)
    dw2 <- sqrt(dt)*rnorm(n)
    x_sim <- rep(x0,n)
    y_sim <- rep(y0,n)
    for (s in 2:n){
      # create anthithetic paths
      x_sim[s] <- x_sim[s-1]-a*x_sim[s-1]*dt+sig*dw1[s]
      y_sim[s] <- y_sim[s-1]-b*y_sim[s-1]*dt+eta*(rho*dw1[s]+sqrt(1-rho^2)*dw2[s])
    }
    sim_i <- c(r0, x_sim[2:n]+y_sim[2:n]+phi_t)
    sim <- append(sim, sim_i)
  }
  out <- matrix(unlist(sim), nSim, n, byrow = T)
  return(out)  
}
##### BOND PRICING #####
# Zero-coupon bond
zcBond <- function(paths, fv){
  dt <- 1/360
  pv <- fv*mean(exp(-apply(paths[,2:ncol(paths)]*dt, 1, sum)))
  return(pv)
}
# Coupon bond
cpnBond <- function(paths, t, c, fv){
  pmtT <- seq(0.5, t, 0.5)*360
  pmtC <- c(rep(c, length(pmtT)-1), c+fv)
  pmtPV <- c()
  for (i in 1:length(pmtT)){
    pmtPV[i] <- zcBond(paths[,1:pmtT[i]],pmtC[i])
  }
  return(sum(pmtPV))
}
# explicit Vasicek zero-coupon bond price
zcBond_vasicek <- function(paths, sig, kappa, r_bar, exerT, t, fv){
  exerT_day <- exerT*360+1
  rt <- paths[,exerT_day]
  B <- 1/kappa*(1-exp(-kappa*(t-exerT)))
  A <- exp((r_bar - sig^2/(2*kappa^2))*(B-(t-exerT)) - sig^2/(4*kappa)*B^2)
  return(fv*mean(A*exp(-B*rt)))
}
##### OPTION PRICING #####
# European call option on ZC bond
zcOption_vasicek <- function(paths, sig, kappa, r_bar, exerT, t, k, fv, type){
  dt <- 1/360
  exerT_day <- exerT*360+1
  rt <- paths[,exerT_day]
  B <- 1/kappa*(1-exp(-kappa*(t-exerT)))
  A <- exp((r_bar - sig^2/(2*kappa^2))*(B-(t-exerT)) - sig^2/(4*kappa)*B^2)
  bond_exerT <- fv*(A*exp(-B*rt))
  exerValue <- switch (type,
    "call" = ifelse(bond_exerT>k, bond_exerT-k, 0),
    "put"  = ifelse(bond_exerT<k, k-bond_exerT, 0)
  )
  option <- mean(exerValue*exp(-apply(paths[,2:exerT_day]*dt, 1, sum)))
  return(option)
}
# European call option on coupon bond
cpnOption <- function(paths, exerT, t, k, fv, type, method){
  dt <- 1/360
  exerT_to_t <- (t-exerT)*360
}