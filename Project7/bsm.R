# Black-schole pricing model for European options
bsCall <- function(s0,t,x,r,sigma){
  # compute d1, d2
  d1 <- (log(s0/x)+(r+sigma^2/2)*t)/(sigma*sqrt(t))
  d2 <- d1-sigma*sqrt(t)
  # output: Black-Sholes price
  c_bs <- s0*pnorm(d1)-x*exp(-r*t)*pnorm(d2) 
  return(c_bs)
}
bsPut <- function(s0,t,x,r,sigma){
  # compute d1, d2
  d1 <- (log(s0/x)+(r+sigma^2/2)*t)/(sigma*sqrt(t))
  d2 <- d1-sigma*sqrt(t)
  # output: Black-Sholes price
  p_bs <- x*exp(-r*t)*pnorm(-d2)-s0*pnorm(-d1) 
  return(p_bs)
}