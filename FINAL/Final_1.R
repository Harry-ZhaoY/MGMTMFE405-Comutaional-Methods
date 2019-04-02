# function
q1_payoff <- function(x){
  s <- 0
  counter <- 0
  while (s <= x) {
    s <- s+runif(1)
    counter <- counter + 1 
  }
  return(counter)
}
# simulation
x <- 1.1
sim <- c()
for(i in 1:10^4){
  sim[i] <- q1_payoff(x)
}
out <- mean(ifelse(sim>4.54, sim-4.54, 0))