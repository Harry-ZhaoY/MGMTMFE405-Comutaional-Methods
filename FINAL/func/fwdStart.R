# function to price a forward-start European put option
fwdPutEu <- function(s0, r, sigma, t_x, t_n, nSim, n){
  dt <- t_n/(n+1)
  # find # steps for t
  fwdStep <- ceiling(t_x/dt)
  paths <- sPaths(s0, r, sigma, t_n, nSim, n)
  put <- exp(-r*t_n)*mean(ifelse(paths[,fwdStep]>paths[,n], paths[,fwdStep]-paths[,n], 0))
  return(put)
}
# function to price a forward-start American put option
fwdPutAm <- function(s0, r, sigma, t_x, t_n, nSim, n, func, term){
  dt <- t_n/(n+1)
  fwdStep <- ceiling(t_x/dt)
  paths <- sPaths(s0, r, sigma, t_n, nSim, n)
  strike <- paths[,fwdStep]
  # use LSMC to find put price at t_x using stock price from t_x to t_n
  put_tx <- lsmc(paths[,fwdStep:n], strike, r, t_n-t_x, nSim, n-fwdStep+1, func, term)
  put <- exp(-r*t_x)*put_tx
  return(put)
}
