fsLookback <- function(s0, x, r, sigma, t, nSim, n, type){
  paths <- sPaths(s0, r, sigma, t, nSim, n)
  if (type == "call"){
    maxS <- apply(paths, 1, max)
    v <- exp(-r*t)*mean(ifelse(maxS>x, maxS-x, 0))
  }
  else if (type == "put"){
    minS <- apply(paths, 1, min)
    v <- exp(-r*t)*mean(ifelse(x>minS, x-minS, 0))
  }
  else {
    stop("Incorrect option type!")
  }
  return(v)
}