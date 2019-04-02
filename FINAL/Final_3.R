##### Final 3 #####
# stock paths
sPaths <- function(s0 = 100, r = 0.05, sigma = 0.35, t = 5, nSim = 1000, n = 1000){
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
# function to find the price of this specified option
thisOption <- function(paths, r = 0.05, k = 100, termT = 5, n = 1000){
  t <- seq(0,termT,termT/n)
  L <- 50*exp(0.138629*t)
  U <- 200-50*exp(0.138629*t)
  taoL <- apply(paths, 1, function(path){return(which(path<=L)[1])})
  taoU <- apply(paths, 1, function(path){return(which(path>=U)[1])})
  taoL[which(is.na(taoL))] <- taoU[which(is.na(taoU))]<- n
  payoff <- c()
  for (i in 1:nrow(paths)){
    if (taoL[i]==taoU[i]){
      payoff[i] <- 0
    } else{
      if(taoL[i]<taoU[i]){
        payoff[i] <- exp(-r*(taoL[i]/n)*termT)*(k-paths[i,taoL[i]])
      } else {
        payoff[i] <- exp(-r*(taoU[i]/n)*termT)*(paths[i,taoU[i]]-k)
      }
    }
  }
  return(list(p = payoff, u = taoU, l = taoL))
}

# main
thisPaths <- sPaths()
out <- thisOption(thisPaths)
# a. price of the option
mean(out$p)
# b. prob
# number of exercise
n <- sum(out$p!=0)
# the conditional prob
p <- sum(out$l<out$u)/n
