# Project 3
##### FUNCTIONS FROM PROJECT 1 #####
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

##### FUNCTIONS #####
## Question 1 ##
# function to simulate x2
q1_x <- function(seed, x0, intv, n, dt){
  # create a matrix of normal random variables
  dWt <- sqrt(dt)*matrix(sample(rnorm_pmm(seed, 10^5), intv[2]/dt*n, T),intv[2]/dt,n) 
  x_sim <- rep(0,n)
  for (i in 1:n){
    x <- rep(x0, intv[2]/dt)
    for (t in 2:(intv[2]/dt)){
      x[t] <- x[t-1] + (1/5-1/2*x[t-1])*dt + 2/3*dWt[t,i]
    }
    x_sim[i] <- x[t]
  }
  # return a vector of simulated x2 values
  return(x_sim)
}
# function to simulate yt
q1_y <- function(seed, y0, intv, n, dt){
  dZt <- sqrt(dt)*matrix(sample(rnorm_pmm(seed, 10^5), intv[2]/dt*n, T),intv[2]/dt,n) 
  y_sim <- rep(0,n)
  seq_t <- seq(intv[1],intv[2],dt)
  for (i in 1:n){
    y <- rep(y0, intv[2]/dt)
    for (t in 2:(intv[2]/dt)){
      y[t] <- y[t-1] + (2/(1+seq_t[t])*y[t-1] + (1+seq_t[t]^3)/3)*dt + (1+seq_t[t]^3)/3*dZt[t,i]
    }
    y_sim[i] <- y[t]
  }
  # return a vector of simulated x2 values
  return(y_sim)
}
# cube root function to handle NaN
cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}
## Question 2 ##
# function to find x3
q2_x <- function(seed1,seed2, x0, intv, n, dt){
  # generate random variables
  dWt <- sqrt(dt)*matrix(rnorm(intv[2]/dt*n),intv[2]/dt,n) 
  dZt <- sqrt(dt)*matrix(rnorm(intv[2]/dt*n),intv[2]/dt,n) 
  # initiate xt
  x_sim <- rep(0,n)
  for (i in 1:n){
    x <- rep(x0, intv[2]/dt)
    for (t in 2:(intv[2]/dt)){
      x[t] <- x[t-1] + 1/4*x[t-1]*dt + 1/3*x[t-1]*dWt[t,i] - 3/4*x[t-1]*dZt[t,i]
    }
    x_sim[i] <- sum(x*dt) 
  }
  # return a vector of simulated x2 values
  return(x_sim)
}
## Question 3 ##
# monte carlo simulation function
callprice <- function(seed, s0, t, x, r, sigma){
  n <- 2^26
  st <- s0*exp((r-sigma^2/2)*t + sigma*sqrt(t)*rnorm(n))  
  c <- exp(-r*t)*mean(ifelse((st-x)>0,st-x,0))
  return(c)
}
callprice_bs <- function(s0,t,x,r,sigma){
  # compute d1, d2
  d1 <- (log(s0/x)+(r+sigma^2/2)*t)/(sigma*sqrt(t))
  d2 <- d1-sigma*sqrt(t)
  # output: Black-Sholes price
  cprice_bs <- s0*pnorm(d1)-x*exp(-r*t)*pnorm(d2) 
  return(cprice_bs)
}
## Question 4 ##
# 2-factor model function
q4 <- function(seed1, seed2){
  dt <- 0.01
  n <- 1000
  # generate random variables
  z1 <- sample(rnorm_pmm(seed1, 10^5), t/dt*n, T)
  z2 <- sample(rnorm_pmm(seed2, 10^5), t/dt*n, T)
  # correlated processes
  dWt1 <- matrix(sqrt(dt)*z1,t/dt,n) 
  dWt2 <- matrix(sqrt(dt)*(rho*z1 + sqrt(1-rho^2)*z2),t/dt,n)
  # simulation
  s_sim <- list(ft = rep(0,n), 
                pt = rep(0,n), 
                re = rep(0,n))
  for (i in 1:n){
    st_ft <- rep(s0, t/dt)
    vt_ft <- rep(v0, t/dt)
    st_pt <- rep(s0, t/dt)
    vt_pt <- rep(v0, t/dt)
    st_re <- rep(s0, t/dt)
    vt_re <- rep(v0, t/dt)
    for (s in 2:(t/dt)){
      # full truncation
      vt_ft[s] <- vt_ft[s-1] + alpha*(beta - ifelse(vt_ft[s-1]>=0, vt_ft[s-1], 0))*dt + sigma*sqrt(ifelse(vt_ft[s-1]>=0, vt_ft[s-1], 0))*dWt2[s,i]
      st_ft[s] <- st_ft[s-1] + r*st_ft[s-1]*dt + sqrt(ifelse(vt_ft[s-1]>=0, vt_ft[s-1], 0))*st_ft[s-1]*dWt1[s,i]
      # partial truncation
      vt_pt[s] <- vt_pt[s-1] + alpha*(beta - vt_pt[s-1])*dt + sigma*sqrt(ifelse(vt_pt[s-1]>=0, vt_pt[s-1], 0))*dWt2[s,i]
      st_pt[s] <- st_pt[s-1] + r*st_pt[s-1]*dt + sqrt(ifelse(vt_pt[s-1]>=0, vt_pt[s-1], 0))*st_pt[s-1]*dWt1[s,i]
      # reflection
      vt_re[s] <- abs(vt_re[s-1]) + alpha*(beta - abs(vt_re[s-1]))*dt + sigma*sqrt(abs(vt_re[s-1]))*dWt2[s,i]
      st_re[s] <- st_re[s-1] + r*st_re[s-1]*dt + sqrt(abs(vt_re[s-1]))*st_re[s-1]*dWt1[s,i]
    }
    s_sim$ft[i] <- st_ft[s]
    s_sim$pt[i] <- st_pt[s]
    s_sim$re[i] <- st_re[s]
  }
  return(s_sim)
}
## Question 5 ##
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
# question 5e estimate integral
q5e <- function(num, base){
  x <- haltonSeq(num, base[1])
  y <- haltonSeq(num, base[2])
  output <- sum(exp(-x*y)*(sin(6*pi*x) + cbrt(cos(2*pi*y)))*(1/num))
  return(output)
}
###### MAIN #####
seed1 <- as.numeric(Sys.time())*10^5
seed2 <- as.numeric(Sys.time())^2
##### Question 1 #####
# run simulations
x2 <- q1_x(seed1, 1, c(0,2), 1000, 0.01)
y2 <- q1_y(seed1, 3/4, c(0,2), 1000, 0.01)
y3 <- q1_y(seed1, 3/4, c(0,3), 1000, 0.01)
# outputs
prob <- sum(y2[y2>5])/n
E1 <- mean(cbrt(x2))
E2 <- mean(y3)
E3 <- mean((x2*y2)[x2>1])

##### Question 2 #####
# outputs
E1 <- mean((1+q2_x(seed1, seed2, 1, c(0,3), 1000, 0.01))^(1/3))
E2 <- mean((1+exp(-0.08*3+1/3*sqrt(3)*rnorm_pmm(seed1, 10^3) + 3/4*sqrt(3)*rnorm_pmm(seed2, 10^5)))^(1/3))

##### Question 3 #####
# (a): compute call option price with monte carlo simulation
# set parameters
s0 <- 20
t <- 0.5
x <- 20
r <- 0.04
sigma <- 0.25
# output: monte carlo price
C1 <- callprice(seed1,s0,t,x,r,sigma)
# (b): compute call option price with black-scholes

C2 <- callprice_bs(s0,t,x,r,sigma)
# (c) greeks
# set initial prices
s0 <- 15:25
# outputs
delta <- rep(0,11)
theta <- rep(0,11)
vega  <- rep(0,11)
rho   <- rep(0,11)
gamma <- rep(0,11)
for (i in 1:11){
  del <- 0.00005
  cprice <- callprice(seed3, s0[i], t, x, r, sigma)
  delta[i] <- (callprice(seed3, s0[i]+del, t, x, r, sigma)-cprice)/del
  theta[i] <- (callprice(seed3, s0[i], t-del, x, r, sigma)-cprice)/del
  vega[i]  <- (callprice(seed3, s0[i], t, x, r, sigma+del)-cprice)/del
  rho[i]   <- (callprice(seed3, s0[i], t, x, r+del, sigma)-cprice)/del
  gamma[i] <- (callprice(seed3, s0[i]+del, t, x, r, sigma)-2*cprice+callprice(seed3, s0[i]-del, t, x, r, sigma))/(del)^2
} 
# plots
plot(delta, type = "l", lwd = 3, col = "dodgerblue4", main = "Delta")
plot(theta, type = "l", lwd = 3, col = "firebrick2", main = "Theta")
plot(vega, type = "l", lwd = 3, col = "skyblue1", main = "Vega")
plot(gamma, type = "l", lwd = 3, col = "orchid2", main = "Gamma")
plot(rho, type = "l", lwd = 3, col = "palegreen3", main = "Rho")

##### Question 4 #####
# set parameters
rho   <- -0.6
r     <- 0.03
s0    <- 48
v0    <- 0.05
sigma <- 0.42
alpha <- 5.8
beta  <- 0.0625
t     <- 0.5 
x     <- 50
# outputs
outQ4 <- q4(seed1, seed2)
C1 <- mean(ifelse(outQ4$ft-x>0, outQ4$ft-x, 0))
C2 <- mean(ifelse(outQ4$pt-x>0, outQ4$pt-x, 0))
C3 <- mean(ifelse(outQ4$re-x>0, outQ4$re-x, 0))

##### Question 5 #####
# (a) generate 100 2-dimension vector of U[0,1]x[0,1]
# uniform number generation
u100x2 <- cbind(runif(seed1, 100), runif(seed2*5, 100))
# (b) halton sequences base(2,7)
h100x2a <- cbind(haltonSeq(100, 2), haltonSeq(100, 7))
# (c) halton sequences base(2,4)
h100x2b <- cbind(haltonSeq(100, 2), haltonSeq(100,4))
# (d) plot the sequences
plot(u100x2[,1],u100x2[,2])
plot(h100x2a[,1],h100x2a[,2])
plot(h100x2b[,1],h100x2b[,2])
# (e) 
q5e(10000, c(2,4))
q5e(10000, c(2,7))
q5e(10000, c(5,7))

