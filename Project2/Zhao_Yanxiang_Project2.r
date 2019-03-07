# Project 2
# FUNCTIONS FROM PROJECT 1
# LGM random number generation
runif <- function(seed, n){
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
rnorm_pmm <- function(seed, n){
  # define the output vector
  output <- c()
  trial <- 1
  while(length(output)<n){
    m <- n - length(output)
    u0 <- runif(seed*trial,m)
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
    trial <-  trial + 1
  }
  return(output)
}

# ANSWERS
# set seed
seed <<- as.numeric(Sys.time())
##### Question 1 #####
q1 <- function(seed){
  n <- 1000
  z1 <- rnorm_pmm(seed, n)
  z2 <- rnorm_pmm(seed*100, n)
  mu <- c(0,0)
  sigma <- matrix(c(1,-0.7, -0.7, 1), 2, 2)
  # generate x~N(0,1) and y~N(0,1) with Cov(x,y) = -0.7
  x <- mu[1] + sigma[1,1]*z1
  y <- mu[2] + sigma[1,2]/sigma[1,1]*z1 + sigma[2,2]*sqrt(1-(sigma[1,2]/(sigma[1,1]*sigma[2,2]))^2)*z2
  # compute correlation
  rho <- (1/(n-1)*sum((x-sum(x)/n)*(y-sum(y)/n)))/
         (sqrt(1/(n-1)*sum((x-mean(x))^2))*sqrt(1/(n-1)*sum((y-mean(y))^2)))
  return(rho)
}
# find the rho
rho <- q1(seed)

##### Question 2 #####
q2 <- function(seed){
  n <- 1000
  z1 <- rnorm_pmm(seed, n)
  z2 <- rnorm_pmm(seed*100, n)
  rho <- 0.6
  # generate x~N(0,1) and y~N(0,1) with Cov(x,y) = -0.7
  x <- z1
  y <- rho*z1 + sqrt(1-rho^2)*z2
  # compute expected value
  output <- mean(apply(cbind(rep(0, n), x^3+sin(y)+x^2*y),1,max))
  return(output)
}
# find the expected value
E <- q2(seed)

##### Question 3 #####
# (a) Estimate expected values by simulation
q3a1 <- function(seed){
  w_5 <- sqrt(5)*rnorm_pmm(seed, 1000)
  output <- mean(w_5^2+sin(w_5))
  return(output)
}
q3a2 <- function(seed, t){
  z1 <- rnorm_pmm(seed, 1000)
  w <- sqrt(t)*z1
  output <- mean(exp(t/2)*cos(w))
  return(output)
}
# outputs
Ea1 <- q3a1(seed)
Ea2 <- q3a2(seed, 0.5)
Ea3 <- q3a2(seed, 3.2)
Ea4 <- q3a2(seed, 6.5)

# (b) Refer to PDF

# (c) variance reduction
q3b1 <- function(seed){
  z1 <- rnorm_pmm(seed, 1000)
  w_5a <- sqrt(5)*z1
  w_5b <- sqrt(5)*(-z1)
  output <- mean(apply(cbind(w_5a^2+sin(w_5a),w_5b^2+sin(w_5b)), 1, mean))
  return(output)
}
q3b2 <- function(seed, t){
  z1 <- rnorm_pmm(seed, 1000)
  w_a <- sqrt(t)*z1
  w_b <- sqrt(t)*(-z1)
  output <- mean(apply(cbind(exp(t/2)*cos(w_a),exp(t/2)*cos(w_b)),1,mean))
  return(output)
}
# outputs: variance reduction
Eb1 <- q3b1(seed)
Eb2 <- q3b2(seed, 0.5)
Eb3 <- q3b2(seed, 3.2)
Eb4 <- q3b2(seed, 6.5)

##### Question 4 #####
# (a) Estimate the price of a European Call option with simulation
q4a <- function(seed){
  n <- 1000
  # define variables as given
  r <- 0.04
  sigma <- 0.2
  s_0 <- 88
  t <- 5
  x <- 100
  # generate random normal variates
  z1 <- rnorm_pmm(seed, n)
  # calculate the payoffs
  s_T <- s_0*exp((r-sigma^2/2)*t+sigma*sqrt(t)*z1)
  payoffs <- apply(cbind(rep(0,n), s_T-x), 1, max)
             # apply the max function each row to find payoffs >= 0
  # take PV of expected payoff
  c <- exp(-r*t)*mean(payoffs)
  return(c)
}
# output
Ca1 <- q4a(seed)

# (b) 1. Compute option price by Black-Sholes formula
# define variables as given
r <- 0.04
sigma <- 0.2
s_0 <- 88
t <- 5
x <- 100
# compute d1, d2
d1 <- (log(s_0/x)+(r+sigma^2/2)*t)/(sigma*sqrt(t))
d2 <- d1-sigma*sqrt(t)
# output: Black-Sholes
Cb1 <- s_0*pnorm(d1)-x*exp(-r*t)*pnorm(d2) 
       #pnorm() is the built-in normal cumulative density function
# (b) 2. use variance reduction
q4b <- function(seed){
  n <- 1000
  # define variables as given
  r <- 0.04
  sigma <- 0.2
  s_0 <- 88
  t <- 5
  x <- 100
  # generate random normal variates
  z1 <- rnorm_pmm(seed, n)
  # calculate the payoffs with antithetic variates
  # (+z1)
  s_T1 <- s_0*exp((r-sigma^2/2)*t+sigma*sqrt(t)*z1)
  payoffs1 <- apply(cbind(rep(0,n), s_T1-x), 1, max)
  c1 <- exp(-r*t)*(payoffs1)
  # (-z1)
  s_T2 <- s_0*exp((r-sigma^2/2)*t+sigma*sqrt(t)*(-z1))
  payoffs2 <- apply(cbind(rep(0,n), s_T2-x), 1, max)
  c2 <- exp(-r*t)*(payoffs2)
  # take the expected value of the mean
  c <- mean(apply(cbind(c1,c2),1,mean))
  return(c)
}
# output
Cb2 <- q4b(seed)

##### Question 5 #####
# (a) find E[Sn] and plot 
q5a <- function(seed, t){
  n <- 1000
  # define variables as given
  r <- 0.04
  sigma <- 0.18
  s_0 <- 88
  # generate random normal variates
  z1 <- rnorm_pmm(seed, n)
  # calculate s_T with antithetic variates
  # (+z1)
  s_T1 <- s_0*exp((r-sigma^2/2)*t+sigma*sqrt(t)*z1)
  # (-z1)
  s_T2 <- s_0*exp((r-sigma^2/2)*t+sigma*sqrt(t)*(-z1))
  # take the expected value of the mean
  s_T <- mean(apply(cbind(s_T1,s_T1),1,mean))
  return(s_T)
}
# generate the path
path <- rep(0,10)
for(t in 1:10){
  path[t] <- q5a(seed*t, t)
}
plot(path, type = "l")

# (b) finer time intervals
N <- 1000
# generate 6 paths
paths <- matrix(rep(rep(0,N),6),N,6)
for (i in 1:6) {
  step <- 1
  # take finer time intervals 
  for(t in seq(1, 10, (10-1)/(N-1))){
    paths[step,i] <- q5a(seed*t*i, t)
    step <- step + 1
  }
}

# (c) plot a and b in one plot
time <- matrix(rep(seq(1, 10, (10-1)/(N-1)),6),N,6)
matplot(time, paths, type = "l")
par(new=TRUE)
plot(path, type = "l", lwd = 4, col = "red", 
     ylim=c(min(paths), max(paths)), xlab = "time", ylab = "paths")

# (d) refer to PDF

##### Problem 6 #####
# (a) Euler method
q6a <- function(x){
  pi <- 4*sqrt(1-x^2)
  return(pi)
}
# Riemann sum
rSum <- 0
loops <- 1000000
for(i in seq(0,1,1/loops)){
  rSum <- rSum + q6a(i)*(1/loops)
}
# output
Ia <- rSum

# (b) simulation
q6b <- function(seed){
  pi <- mean(4*sqrt(1-runif(seed, 1000)^2))
  return(pi)
}
Ib <- q6b(seed)

# (c) importance sampling method
# model by Wei Cai, from the notes
x <- runif(seed, 1000)
a <- 0.74
g_x <- 4*sqrt(1-x^2)
t_x <- (1-0.74*x^2)/(1-a/3)
# output
Ic <- mean(g_x/t_x)

