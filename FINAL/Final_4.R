##### Final 4 #####
# function to find prices of a zero coupon bond
zcBond <- function(paths, fv = 10000, t = 0.5, n = 100){
  dt <- t/n
  # if paths is a matrix
  if (!is.null(ncol(paths))){
    pv <- fv*(exp(-apply(paths[,2:ncol(paths)]*dt, 1, sum)))
  } else{
    # deal with a single paths
    pv <- fv*(exp(sum(paths[2:length(paths)]*dt)))
  }
  return(pv)
}
# function to simulate interest rate paths
rPaths <- function(r0 = 0.05, alpha = 0.36, beta = -5.86, sig = 0.36, gamma = 2, t = 0.5, nSim = 1000, n = 100){
  dt <- t/n 
  sim <- list()
  for (i in 1:(nSim/2)){
    # generate Brownian motion for the path
    dw <- sqrt(dt)*rnorm(n)
    # simulate the rate paths
    r_sim1 <- rep(r0,n)
    r_sim2 <- rep(r0,n)
    for (s in 2:n){
      # create anthithetic paths
      r_sim1[s] <- r_sim1[s-1] + (alpha+beta*r_sim1[s-1])*dt + sig*r_sim1[s-1]^gamma*dw[s-1]
      r_sim2[s] <- r_sim2[s-1] + (alpha+beta*r_sim2[s-1])*dt - sig*r_sim2[s-1]^gamma*dw[s-1]
    }
    sim_i <- list(r_sim1, r_sim2)
    sim <- append(sim, sim_i)
  }
  out <- matrix(unlist(sim), nSim, n, byrow = T)
  return(out)  
}
# price european put on zero coupon bonds
zcPut <- function(paths_before, k = 9800, exerT = 0.5, termT = 1, mSim = 100){
  option <- c()
  for (i in 1:ncol(paths_before)){
    paths_after <- rPaths()
    bond_exerT <- zcBond(paths_after)
    payoff <- ifelse(bond_exerT<k, k-bond_exerT, 0)
    option[i] <- mean(zcBond(paths_before[i,], payoff))
  }
  return(option)
}
# main
paths <- rPaths()
mean(zcPut(paths))
