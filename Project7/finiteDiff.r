# Zhao_Yanxiang_Projcet7
########################### FUNCTIONS ########################### 
# Explicit Finite-Difference Method 
efd <- function(type, euro, s0, k, r, sig, t, dt, dsx, log){
  # if change in log price
  if (log) {
    n <- floor(log(s0)/dsx)
    # defind M and N
    M <- t/dt+1
    N <- 2*n+1
    # define dx
    dx <- dsx
    # find terminal conditions
    sT <- exp(c(log(s0)+dx*(n:-n)))
    # calculate pu, pm, pd
    pu <- dt*(sig^2/(2*dx^2) + (r-sig^2/2)/(2*dx))
    pm <- 1 - dt*sig^2/dx^2 - r*dt
    pd <- dt*(sig^2/(2*dx^2) - (r-sig^2/2)/(2*dx))
    # construct matrix A
    A <- diag(pm, N-2, N-2)
    diag(A[-1,]) <- pu
    diag(A[,-1]) <- pd
    # set offset constans
    b <- c(pu, pd)
  } else { # if change in price
    # defind M
    M <- t/dt+1
    # define ds
    ds <- dsx
    # find terminal conditions
    sT <- seq(2*s0, 0, -ds)
    N <- length(sT)
    # calculate pu, pm, pd
    j <- (N-2):1
    pu <- 0.5*dt*(sig^2*j^2 + r*j)
    pm <- 1-dt*(sig^2*j^2+r)
    pd <- 0.5*dt*(sig^2*j^2 - r*j)
    # set A
    A <- diag(pm)
    diag(A[-1,]) <- pu[2:length(pu)]
    diag(A[,-1]) <- pd[1:length(pd)-1]
    # set offset constans
    b <- c(pu[1], pd[length(pd)])
  }
  # payoff matrix
  payoff <- diag(0, N, M)
  if (type == "call"){
    payoff[,M] <- ifelse(sT-k>0, sT-k, 0)
    payoff[1,] <- (max(sT)-k)*exp(-r*seq(t,0,-dt))
    payoff[N,] <- 0
  } 
  else if (type == "put"){
    payoff[,M] <- ifelse(k-sT>0, k-sT, 0)
    payoff[1,] <- 0
    payoff[N,] <- (k-min(sT))*exp(-r*seq(t,0,-dt))
  } 
  else{
    stop("Incorrect option type!")
  }
  # for loop through M
  for (i in (M-1):1){
    # if Europian option
    if (euro){
      payoff[2:(N-1),i] <- A%*%payoff[2:(N-1),i+1]
      payoff[2,i] <- payoff[2,i]+b[1]*payoff[1,i+1]
      payoff[N-1,i] <- payoff[N-1,i]+b[2]*payoff[N,i+1]
    } else { # if American option
      cv <- payoff[,i]
      cv[2:(N-1)] <- A%*%payoff[2:(N-1),i+1]
      cv[2] <- payoff[2,i]+b[1]*payoff[1,i+1]
      cv[N-1] <- payoff[N-1,i]+b[2]*payoff[N,i+1]
      payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
    }
  }
  return(payoff[floor(N/2)+1,1])
}
# Implicit Finite-Difference Method
ifd <- function(type, euro, s0, k, r, sig, t, dt, dsx, log){
  # if change in log price
  if (log) {
    n <- floor(log(s0)/dsx)
    # defind M and N
    M <- t/dt+1
    N <- 2*n+1
    # define dx
    dx <- dsx
    # find terminal conditions
    sT <- exp(c(log(s0)+dx*(n:-n)))
    # calculate pu, pm, pd
    pu = -0.5*dt*(sig^2/dx^2 + (r-sig^2/2)/dx)
    pm = 1 + dt*sig^2/dx^2 + r*dt
    pd = -0.5*dt*(sig^2/dx^2 - (r-sig^2/2)/dx)
    # construct matrix A
    A <- diag(pm, N-2, N-2)
    diag(A[-1,]) <- pu
    diag(A[,-1]) <- pd
    A <- rbind(c(1, -1, rep(0,N-2)), 
               cbind(c(pu, rep(0,N-3)),A,c(rep(0,N-3), pd)), 
               c(rep(0,N-2), 1, -1))
  } else { # if change in price
    # defind M
    M <- t/dt+1
    # define ds
    ds <- dsx
    # find terminal conditions
    sT <- seq(2*s0, 0, -ds)
    N <- length(sT)
    # calculate pu, pm, pd
    j <- (N-2):1
    pu <- -0.5*dt*(sig^2*j^2 + r*j)
    pm <- 1+dt*(sig^2*j^2+r)
    pd <- 0.5*dt*(-sig^2*j^2 + r*j)
    # set A
    A <- diag(pm)
    diag(A[-1,]) <- pu[2:length(pu)]
    diag(A[,-1]) <- pd[1:length(pd)-1]
    # set offset constans
    b <- c(-pu[1], -pd[length(pd)])
  }
  # payoff matrix
  payoff <- diag(0, N, M)
  if (type == "call"){
    payoff[,M] <- ifelse(sT-k>0, sT-k, 0)
    payoff[1,] <- (max(sT)-k)*exp(-r*seq(t,0,-dt))
    payoff[N,] <- 0
  } 
  else if (type == "put"){
    payoff[,M] <- ifelse(k-sT>0, k-sT, 0)
    payoff[1,] <- 0
    payoff[N,] <- (k-min(sT))*exp(-r*seq(t,0,-dt))
  } 
  else{
    stop("Incorrect option type!")
  }
  # find A inverse
  Ainv <- solve(A)
  # find option price
  if (log) {
    for (i in (M-1):1){
      # construct matrix B
      B <- payoff[,i+1]
      # if Europian option
      if (euro){
        payoff[,i] <- Ainv%*%B
      } else { # if American option
        cv <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      }
    }
  } else {
    # for loop through M
    for (i in (M-1):1){
      # construct matrix B
      B <- payoff[2:(N-1),i+1]
      B[1] <- B[1] + payoff[1,i+1]*b[1]
      B[length(N)] <- B[length(N)] + payoff[N,i+1]*b[2]
      # if Europian option
      if (euro){
        payoff[2:(N-1),i] <- Ainv%*%B
      } else { # if American option
        cv <- payoff[,i]
        cv[2:(N-1)] <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      }
    }
  }
  return(payoff[floor(N/2)+1,1])
}
# Crank-Nicolson Finite-Difference Method
cnfd <- function(type, euro, s0, k, r, sig, t, dt, dsx, log){
  # if change in log price
  if (log) {
    n <- floor(log(s0)/dsx)
    # defind M and N
    M <- t/dt+1
    N <- 2*n+1
    # define dx
    dx <- dsx
    # find terminal conditions
    sT <- exp(c(log(s0)+dx*(n:-n)))
    # calculate pu, pm, pd
    pu = -0.25*dt*(sig^2/dx^2 + (r-sig^2/2)/dx)
    pm = 1 + dt*sig^2/(2*dx^2) + r*dt/2
    pd = -0.25*dt*(sig^2/dx^2 - (r-sig^2/2)/dx)
    # construct matrix A
    A <- diag(pm, N-2, N-2)
    diag(A[-1,]) <- pu
    diag(A[,-1]) <- pd
    A <- rbind(c(1, -1, rep(0,N-2)), 
               cbind(c(pu, rep(0,N-3)),A,c(rep(0,N-3), pd)), 
               c(rep(0,N-2), 1, -1))
  } else { # if change in price
    # defind M
    M <- t/dt+1
    # define ds
    ds <- dsx
    # find terminal conditions
    sT <- seq(2*s0, 0, -ds)
    N <- length(sT)
    # calculate pu, pm, pd
    j <- (N-2):1
    pu <- 0.25*dt*(sig^2*j^2 + r*j)
    pm <- -0.5*dt*(sig^2*j^2+r)
    pd <- 0.25*dt*(sig^2*j^2 - r*j)
    # set C
    A <- diag(1-pm)
    diag(A[-1,]) <- -pu[2:length(pu)]
    diag(A[,-1]) <- -pd[1:length(pd)-1]
    # set D
    D <- diag(1+pm)
    diag(D[-1,]) <- pu[2:length(pu)]
    diag(D[,-1]) <- pd[1:length(pd)-1]
    # set offset constans
    b <- c(pu[1], pd[length(pd)])
  }
  # payoff matrix
  payoff <- diag(0, N, M)
  if (type == "call"){
    payoff[,M] <- ifelse(sT-k>0, sT-k, 0)
    payoff[1,] <- (max(sT)-k)*exp(-r*seq(t,0,-dt))
    payoff[N,] <- 0
  } 
  else if (type == "put"){
    payoff[,M] <- ifelse(k-sT>0, k-sT, 0)
    payoff[1,] <- 0
    payoff[N,] <- (k-min(sT))*exp(-r*seq(t,0,-dt))
  } 
  else{
    stop("Incorrect option type!")
  }
  # find A inverse
  Ainv <- solve(A)
  # find option price
  if (log) {
    for (i in (M-1):1){
      # construct matrix B
      B <- payoff[,i+1]
      B[2:(N-1)] <- -pu*payoff[1:(N-2),i+1] - (pm-2)*payoff[2:(N-1),i+1] - pd*payoff[3:N,i+1] 
      # if Europian option
      if (euro){
        payoff[,i] <- Ainv%*%B
      } else { # if American option
        cv <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      }
    }
  } else {
    # for loop through M
    for (i in (M-1):1){
      # construct matrix B
      B <- D%*%payoff[2:(N-1),i+1]
      B[1] <- B[1] + payoff[1,i+1]*b[1] + payoff[1,i]*b[1]
      B[length(N)] <- B[length(N)] + payoff[N,i+1]*b[2] + payoff[N,i+1]*b[2]
      # if Europian option
      if (euro){
        payoff[2:(N-1),i] <- Ainv%*%B
      } else { # if American option
        cv <- payoff[,i]
        cv[2:(N-1)] <- Ainv%*%B
        payoff[,i] <- apply(cbind(payoff[,M], cv),1,max)
      }
    }
  }
  return(payoff[floor(N/2)+1,1])
}
