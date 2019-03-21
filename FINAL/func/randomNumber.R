runif_lmg <- function(seed, n){
  # initiate the random number series
  x <- c(seed)
  # LGM parameters
  m <- 2^31-1
  a <- 7^5
  b <- 0
  # random number generation
  for(t in 1:n){
    x <- c(x, (a*x[t]+b) %% m)
  }
  # drop the seed
  x <- x[-1]
  # scale to U[0,1]
  return(x/m)
}

# random discrete distribution function
rdisc <- function(rn, val, p){
  # args: rn is a unified randomly distributed series U[0,1],
  #       val is a list of values,
  #       p is a list of probailities, sum to 1
  stopifnot(sum(p)==1) # check if probabilities sum to 1
  stopifnot(length(val)==length(p)) # check if value and prob match
  # generate 0 series
  output <- rep(0, times = length(rn))
  # initial lower-bound
  lb <- 0
  for (i in 1:length(p)){
    # rewrite upper-bound
    ub <- sum(p[1:i])
    # assign discrete value
    output[rn<=ub & rn>lb] <- val[i]
    # rewrite lower-bound
    lb <- ub
  }
  return(output)
}
# random exponentially distributed numbers
rexp <- function(rn, lambda){
  output <- -lambda*log(1-rn)
  return(output)
}
# random normal dist
# Box-Muller method
rnorm_bmm <- function(n){
  m <- n + n %% 2 # to overwirte odd number
  u0 <- runif(m)
  # split into 2 vectors
  u <- split(u0, rep(1:2, m/2))
  # generate N(0,1)
  z1 <- sqrt(-2*log(u[[1]]))*cos(2*pi*u[[2]])
  z2 <- sqrt(-2*log(u[[1]]))*sin(2*pi*u[[2]])
  output <- c(z1, z2)[1:n]
  return(output)
}
# Polar-Marsaglia
rnorm_pmm <- function(n){
  # define the output vector
  output <- c()
  while(length(output)<n){
    m <- n - length(output)
    u0 <- runif(m)
    u <- split(u0, rep(1:2, m/2))
    # Define v1, v2, w
    vw <- as.data.frame(matrix(c(2*u[[1]]-1, 2*u[[2]]-1), m/2, 2))
    vw <- cbind(vw, vw[,1]^2 + vw[,2]^2)
    # select w <= 1
    vw <- vw[vw[,3]<=1,]
    # generate N(0,1)
    z1 <- vw[,1]*sqrt((-2*log(vw[,3]))/vw[,3])
    z2 <- vw[,2]*sqrt((-2*log(vw[,3]))/vw[,3])
    
    output <- c(output, z1, z2)
  }
  return(output)
}