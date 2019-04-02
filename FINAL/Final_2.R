# function for vt and st paths
vsPaths <- function(rho, v0 = 0.06, alpha = 0.45, beta = -5.105, gamma = 0.25, s0 = 20, r = 0.05, t = 2, nSim = 1000, n = 100){
  dt <- t/n
  s_list <- list()
  v_list <- list()
  for (i in 1:nSim){
    dz1 <- rnorm(n)
    dz2 <- rnorm(n)
    # find correlated Brownian motions
    db <- sqrt(dt)*dz1
    dw <- sqrt(dt)*(rho*dz1+sqrt((1-rho^2))*dz2)
    # generate paths
    s_sim <- rep(s0,n+1)
    v_sim <- rep(v0,n+1)
    for (s in 2:(n+1)){
      # vt path
      v_sim[s] <- v_sim[s-1]+(alpha+beta*max(v_sim[s-1],0))*dt + gamma*sqrt(max(v_sim[s-1],0))*db[s-1]
      # stock price path
      s_sim[s] <- s_sim[s-1]*(1+r*dt + sqrt(max(v_sim[s-1],0))*(dw[s-1]))
    }
    s_list <- append(s_list, s_sim)
    v_list <- append(v_list, v_sim)
  }
  s_mat <- matrix(unlist(s_list), nSim, n+1, byrow = T)
  v_mat <- matrix(unlist(v_list), nSim, n+1, byrow = T)
  return(list(st = s_mat, vt = v_mat)) 
}
# function for floating strike asian call
fsAsianCall <- function(paths){
  A <- apply(paths,1,mean)
  ST <- paths[,ncol(paths)]
  payoff <- ifelse(ST>A, ST-A, 0)
  return(mean(payoff))
}

# main
Rho <- c(-0.75,0,0.75)
out <- c()
for (rho in Rho){
  paths <- vsPaths(rho)
  out <- c(out, fsAsianCall(paths$st))
}