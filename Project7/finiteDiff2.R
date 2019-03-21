# dx
explicit_fd = function(s0, sigma, r, dt, time, strike, delta_x_type){
  
  if(delta_x_type == 1){
    dx = sigma*sqrt(dt)
  }
  if(delta_x_type == 2){
    dx = sigma*sqrt(3*dt)
  }
  if(delta_x_type == 3){
    dx = sigma*sqrt(4*dt)
  }
  
  pu = dt*(sigma^2/(2*dx^2) + (r-sigma^2/2)/(2*dx))
  pm = 1 - dt*sigma^2/dx^2 - r*dt
  pd = dt*(sigma^2/(2*dx^2) - (r-sigma^2/2)/(2*dx))
  # number of nodes in th end, col_node stands for nodes from up to down
  # row_node stands for nodes from left to right
  col_node = 2*(time/dt)+1
  row_node = (time/dt)+1
  # following are for generating the A matrix
  A1 = A3 = rep(0, col_node)
  A1[1:3] = c(pu, pm, pd)
  A3[(col_node-2):col_node] = c(pu, pm, pd)
  A2 = A1
  for(i in 1:(col_node-3)){
    
    temp = rep(0, col_node)
    temp[(i+1):(i+3)] = c(pu, pm, pd)
    A2 = rbind(A2, temp)
  }
  A = rbind(A1,A2,A3)
  # uniformly distributed x
  x_seq = c(sapply((col_node%/%2):1, function(x) log(s0)+dx*x), log(s0),
            sapply(1:(col_node%/%2), function(x) log(s0)-dx*x))        
  # stock prices
  stocks = exp(x_seq)
  
  # B matrix
  B = rep(0, col_node)
  B[col_node] = stocks[col_node-1] - stocks[col_node]
  # Fi matrix
  Fi = sapply(stocks, function(x) max(0, strike-x))
  
  
  
  # following for calculating option payoff until reaching t0
  for(i in (row_node):1){
    
    Fi = A%*%Fi + B
  }
  return(Fi[row_node])
}

implicit_fd = function(s0, sigma, r, dt, time, strike, delta_x_type){
  
  if(delta_x_type == 1){
    dx = sigma*sqrt(dt)
  }
  if(delta_x_type == 2){
    dx = sigma*sqrt(3*dt)
  }
  if(delta_x_type == 3){
    dx = sigma*sqrt(4*dt)
  }
  
  pu = -0.5*dt*(sigma^2/dx^2 + (r-sigma^2/2)/dx)
  pm = 1 + dt*sigma^2/dx^2 + r*dt
  pd = -0.5*dt*(sigma^2/dx^2 - (r-sigma^2/2)/dx)
  # number of nodes in th end, col_node stands for nodes from up to down
  # row_node stands for nodes from left to right
  col_node = 2*(time/dt)+1
  row_node = (time/dt)+1
  # following are for generating the A matrix
  A1 = A3 = rep(0, col_node)
  A1[1:2] = c(1, -1)
  A3[(col_node-1):col_node] = c(1, -1)
  A2 = rep(0, col_node)
  A2[1:3] = c(pu, pm, pd)
  for(i in 1:(col_node-3)){
    
    temp = rep(0, col_node)
    temp[(i+1):(i+3)] = c(pu, pm, pd)
    A2 = rbind(A2, temp)
  }
  A = rbind(A1,A2,A3)
  
  # following for populating x uniformally, and then take log to get the price
  x_seq = c(sapply((col_node%/%2):1, function(x) log(s0)+dx*x), log(s0),
            sapply(1:(col_node%/%2), function(x) log(s0)-dx*x))        
  # stock prices
  stocks = exp(x_seq)
  
  # payoff in the end column
  Fi = sapply(stocks, function(x) max(0, strike-x))
  
  # B matrix
  B = Fi
  # S_i,N -1 - S_i,N
  B[col_node] = stocks[col_node-1] - stocks[col_node]
  B[1] = 0
  
  # following for calculating option payoff
  # Fi for put
  invA = solve(A)
  for(i in (row_node):1){
    
    Fi = invA%*%B
    B[2:(col_node-1)] = Fi[2:(col_node-1)]
  }
  return(Fi[row_node])
}

crank_fd = function(s0, sigma, r, dt, time, strike, delta_x_type){
  
  if(delta_x_type == 1){
    dx = sigma*sqrt(dt)
  }
  if(delta_x_type == 2){
    dx = sigma*sqrt(3*dt)
  }
  if(delta_x_type == 3){
    dx = sigma*sqrt(4*dt)
  }
  
  pu = -0.25*dt*(sigma^2/dx^2 + (r-sigma^2/2)/dx)
  pm = 1 + dt*sigma^2/(2*dx^2) + r*dt/2
  pd = -0.25*dt*(sigma^2/dx^2 - (r-sigma^2/2)/dx)
  # number of nodes in th end, col_node stands for nodes from up to down
  # row_node stands for nodes from left to right
  col_node = 2*(time/dt)+1
  row_node = (time/dt)+1
  # following are for generating the A matrix
  A1 = A3 = rep(0, col_node)
  A1[1:2] = c(1, -1)
  A3[(col_node-1):col_node] = c(1, -1)
  A2 = rep(0, col_node)
  A2[1:3] = c(pu, pm, pd)
  for(i in 1:(col_node-3)){
    
    temp = rep(0, col_node)
    temp[(i+1):(i+3)] = c(pu, pm, pd)
    A2 = rbind(A2, temp)
  }
  A = rbind(A1,A2,A3)
  
  # following for populating x uniformally, and then take log to get the price
  x_seq = c(sapply((col_node%/%2):1, function(x) log(s0)+dx*x), log(s0),
            sapply(1:(col_node%/%2), function(x) log(s0)-dx*x))        
  # stock prices
  stocks = exp(x_seq)
  
  # payoff in the end column
  Fi = sapply(stocks, function(x) max(0, strike-x))
  payoff = Fi
  
  # B matrix
  B = Fi
  # S_i,N -1 - S_i,N
  B[col_node] = stocks[col_node-1] - stocks[col_node]
  B[1] = 0
  B[2:(col_node-1)] = -pu*payoff[1:(col_node-2)] - (pm-2)*payoff[2:(col_node-1)] - pd*payoff[3:col_node]
  
  # following for calculating option payoff
  # Fi for put
  invA = solve(A)
  for(i in (row_node):1){
    
    Fi = invA%*%B
    B[2:(col_node-1)] = -pu*Fi[1:(col_node-2)] - (pm-2)*Fi[2:(col_node-1)] - pd*Fi[3:col_node]
  }
  return(Fi[row_node])
}

# ds
# explicit method using ds
explicit_fd_amer = function(s0, sigma, r, dt, time, strike, ds, type){
  
  # period
  n = (time/dt)+1
  # stock price
  stocks = seq(2*s0, 0, -ds)
  m = length(stocks)
  
  # making matrix B, and End payoff for Fi according to different type
  B = rep(0, m)
  Fi = rep(0, m)
  payoff = rep(0, m)
  
  if(type=='call'){
    B[1] = stocks[1]-stocks[2]
    Fi = payoff = sapply(stocks, function(x) max(0, x-strike))
  }
  if(type=='put'){
    B[m] = stocks[m-1] - stocks[m]
    # Fi matrix, the last column payoff
    Fi = payoff = sapply(stocks, function(x) max(0, strike-x))
  }
  
  # following for making matrix A
  pu_end = dt*(r*(m-1)/2 + sigma^2*(m-1)^2/2)
  pm_end = 1 - dt*(sigma^2*(m-1)^2 + r)
  pd_end = dt*(-r*(m-1)/2 + sigma^2*(m-1)^2/2)
  pu_1 = dt*(r*(1)/2 + sigma^2*(1)^2/2)
  pm_1 = 1 - dt*(sigma^2*(1)^2 + r)
  pd_1 = dt*(-r*(1)/2 + sigma^2*(1)^2/2)
  
  A1 = A3 = rep(0, m)
  A1[1:3] = c(pu_end, pm_end, pd_end)
  A3[(m-2):m] = c(pu_1, pm_1, pd_1)
  A2 = A1
  for(i in 2:(m-1)){
    
    pu = dt*(r*(m-i)/2 + sigma^2*(m-i)^2/2)
    pm = 1 - dt*(sigma^2*(m-i)^2 + r)
    pd = dt*(-r*(m-i)/2 + sigma^2*(m-i)^2/2)
    temp = rep(0, m)
    temp[(i-1):(i+1)] = c(pu, pm, pd)
    A2 = rbind(A2, temp)
  }
  A = rbind(A2,A3)
  
  
  # following for calculating option payoff until reaching t0
  for(i in (n):1){
    
    Fi = A%*%Fi + B
    Fi = pmax(Fi, payoff)
  }
  return(Fi[m%/%2+1])
}

# implicit method using ds
implicit_fd_amer = function(s0, sigma, r, dt, time, strike, ds, type){
  
  # period
  n = (time/dt)+1
  # stock price
  stocks = seq(2*s0, 0, -ds)
  m = length(stocks)
  
  # making matrix B, and End payoff for Fi according to different type
  B = rep(0, m)
  Fi = rep(0, m)
  payoff = rep(0, m)
  
  if(type=='call'){
    B = payoff = sapply(stocks, function(x) max(0, x-strike))
    B[1] = stocks[1]-stocks[2]
    B[m] = 0
  }
  if(type=='put'){
    B = payoff = sapply(stocks, function(x) max(0, strike-x))
    B[1] = 0
    B[m] = stocks[m-1] - stocks[m]
  }
  
  # following for making matrix A
  A1 = A3 = rep(0, m)
  A1[1:2] = c(1, -1)
  A3[(m-1):m] = c(1, -1)
  A2 = A1
  for(i in 2:(m-1)){
    
    s = stocks[i]
    pu = -0.5*dt*((sigma*s)^2/ds^2 + (r*s)/ds)
    pm = 1 + dt*(sigma*s)^2/ds^2 + r*dt
    pd = -0.5*dt*((sigma*s)^2/ds^2 - (r*s)/ds)
    
    temp = rep(0, m)
    temp[(i-1):(i+1)] = c(pu, pm, pd)
    A2 = rbind(A2, temp)
  }
  A = rbind(A2,A3)
  
  invA = solve(A)
  # following for calculating option payoff until reaching t0
  for(i in (n):1){
    
    Fi = invA%*%B
    Fi = pmax(Fi, payoff)
    B[2:(m-1)] = Fi[2:(m-1)]
  }
  return(Fi[m%/%2+1])
}

# Crank-Nicolson method using ds
crank_fd_amer = function(s0, sigma, r, dt, time, strike, ds, type){
  
  # period
  n = (time/dt)+1
  # stock price
  stocks = seq(2*s0, 0, -ds)
  m = length(stocks)
  
  # making matrix B, and End payoff for Fi according to different type
  B = rep(0, m)
  Fi = rep(0, m)
  payoff = rep(0, m)
  
  if(type=='call'){
    B = payoff = sapply(stocks, function(x) max(0, x-strike))
    B[1] = stocks[2]-stocks[1]
    B[m] = 0
  }
  if(type=='put'){
    B = payoff = sapply(stocks, function(x) max(0, strike-x))
    B[1] = 0
    B[m] = stocks[m-1] - stocks[m]
  }
  
  # following for making matrix A
  A1 = A3 = rep(0, m)
  A1[1:2] = c(1, -1)
  A3[(m-1):m] = c(1, -1)
  A2 = A1
  pu1 = pm1 = pd1 = vector()
  for(i in 2:(m-1)){
    
    s = stocks[i]
    pu = -0.25*dt*((sigma*s)^2/ds^2 + (r*s)/ds)
    pm = 1 + dt*(sigma*s)^2/(2*ds^2) + r*dt/2
    pd = -0.25*dt*((sigma*s)^2/ds^2 - (r*s)/ds)
    # populating the pu, pm, pd matrix
    pu1[i] = pu
    pm1[i] = pm
    pd1[i] = pd
    
    temp = rep(0, m)
    temp[(i-1):(i+1)] = c(pu, pm, pd)
    A2 = rbind(A2, temp)
  }
  A = rbind(A2,A3)
  # filling the middle part of matrix B
  B[2:(m-1)] = -pu1[2:(m-1)]*payoff[1:(m-2)] - (pm1[2:(m-1)]-2)*payoff[2:(m-1)] - pd1[2:(m-1)]*payoff[3:m]
  # following for calculating option payoff until reaching t0
  invA = solve(A)
  for(i in (n):1){
    
    Fi = invA%*%B
    Fi = pmax(Fi, payoff)
    B[2:(m-1)] = -pu1[2:(m-1)]*Fi[1:(m-2)] - (pm1[2:(m-1)]-2)*Fi[2:(m-1)] - pd1[2:(m-1)]*Fi[3:m]  
    
  }
  
  
  return(Fi[m%/%2+1])
}
