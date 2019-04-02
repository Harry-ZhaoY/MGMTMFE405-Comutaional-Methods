##### Final 4 #####
# stock and exchange rate paths generator
sExPaths <- function(s0=6000, e0=0.0096, r = 0.05, q = 0, rf = 0.04, sig1 = 0.1, sig2 = 0.15, gamma = -0.04, lambda = 1.5, t = 1, nSim = 1000, n = 100){
  dt <- t/n
  s_list <- list()
  e_list <- list()
  rho = -0.25*dt
  for (i in 1:nSim){
    dz1 <- rnorm(n)
    dz2 <- rnorm(n)
    # find correlated Brownian Motions
    dw <- sqrt(dt)*dz1
    db <- sqrt(dt)*(rho*dz1+sqrt((1-rho^2))*dz2)
    # find Possion jumps
    dj <- rpois(n, lambda*dt)
    # paths generateion
    s_sim <- rep(s0,n+1)
    e_sim <- rep(e0,n+1)
    for (s in 2:(n+1)){
      # stock path
      s_sim[s] <- s_sim[s-1]*(1 + (r - q)*dt + sig1*(dw[s-1]) + gamma*dj[s-1])
      # exchange rate paths
      e_sim[s] <- e_sim[s-1]*(1 + (r - rf)*dt + sig2*db[s-1])
    }
    s_list <- append(s_list, s_sim)
    e_list <- append(e_list, e_sim)
  }
  s_mat <- matrix(unlist(s_list), nSim, n+1, byrow = T)
  e_mat <- matrix(unlist(e_list), nSim, n+1, byrow = T)
  return(list(St = s_mat, Et = e_mat)) 
}
# function to price the call option on stock and exchange rate
sExCall <- function(twoPaths, k = 60, r = 0.05, t = 1){
  ST <- twoPaths$St[, ncol(twoPaths$St)]
  ET <- twoPaths$Et[, ncol(twoPaths$Et)]
  payoffs <- ifelse(ST*ET>k, ST*ET-k, 0)
  return(exp(-r*t)*mean(payoffs))
}

# main
# find the paths
paths <- sExPaths()
# find the call price
sExCall(paths)
