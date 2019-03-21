# stock paths
sPaths <- function(s0, r, sigma, t, nSim, n){
  dt <- t/n
  s_sim <- matrix(rep(s0, nSim), ncol = 1)
  for (s in 2:(n+1)){
    dw <- rnorm(nSim)*sqrt(dt)
    # use partial truncation
    newCol <- s_sim[,s-1]*(1 + r*dt + sigma*dw)
    s_sim <- cbind(s_sim, newCol)
  }
  return(s_sim)  
}
# function to simulate callateral values
vPaths <- function(v0, mu, sigma, gamma, t, lambda1, nSim, n){
  dt <- t/n
  sim <- list()
  for (i in 1:(nSim/2)){
    dw <- sqrt(dt)*rnorm(n)
    dj <- rpois(n, lambda1*dt)
    v_sim1 <- rep(v0,n)
    v_sim2 <- rep(v0,n)
    for (s in 2:n){
      # create anthithetic paths
      v_sim1[s] <- v_sim1[s-1]*(1 + mu*dt + sigma*dw[s] + gamma*dj[s])
      v_sim2[s] <- v_sim2[s-1]*(1 + mu*dt - sigma*dw[s] + gamma*dj[s])
    }
    sim_i <- list(v_sim1, v_sim2)
    sim <- append(sim, sim_i)
  }
  out <- matrix(unlist(sim), nSim, n, byrow = T)
  return(out) 
}
# Vasicek model
vasicekPath <- function(r0, sig, kappa, r_bar, t, nSim, days = 360){
  dt <- 1/days # assume 30/360 convension
  n <-  t/dt+1
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
cirPath <- function(r0, sig, kappa, r_bar, t, nSim, days = 360){
  dt <- 1/days # assume 30/360 convension
  n <- t/dt+1
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
g2ppPath <- function(x0, y0, r0, phi_t, rho, a, b, sig, eta, t, nSim, days = 360){
  dt <- 1/days # assume 30/360 convension
  n <- t/dt+1
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
