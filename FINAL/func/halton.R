# Halton's Low-Discrepancy Sequences
haltonSeq <- function(num, base){
  hseq <- rep(0, num)
  numBits <- 1+ceiling(log(num)/log(base))
  vetBase <- base^(-(1:numBits))
  workVet <- rep(0, numBits)
  for(i in 1:num){
    j <- 1
    ok <- 0
    while (ok==0) {
      workVet[j] <- workVet[j]+1
      if (workVet[j]<base){
        ok <-  1
      } else {
        workVet[j] <- 0
        j <- j+1
      }
    }
    hseq[i] <- sum(workVet*vetBase)
  }
  return(hseq)
}
# halton's call option price
call_halton <- function(s0, x, rf, sigma, t, n, b1, b2){
  # get Halton sequences
  h1 <- haltonSeq(n/2, b1)
  h2 <- haltonSeq(n/2, b2)
  # get random normal variables
  z1 <- sqrt(-2*log(h1))*cos(2*pi*h2)
  z2 <- sqrt(-2*log(h1))*sin(2*pi*h2)
  wt <- sqrt(t)*c(rbind(z1, z2)) # c(rbind(z1, z2)) mix z1 and z2
  st <- s0*exp((rf-sigma^2/2)*t + sigma*wt)
  call <- exp(-rf*t)*mean(ifelse((st-x)>0,st-x,0))
  return(call)
}