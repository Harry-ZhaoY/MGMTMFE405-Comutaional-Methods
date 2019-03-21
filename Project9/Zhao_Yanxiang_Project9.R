# Zhao_Yanxiang_Project6
setwd("C:/Users/harry/OneDrive/Documents/GitHub/MGMTMFE405-Comutaional-Methods/Project9/")
source('abs.R')
##### Question 1 #####
# a)
source("promptMbs.R")

# b) for different kappa
kappa <- seq(0.3, 0.9, 0.1)
q1b <- c()
for (k in kappa) {
  value <- findMBS(kappa=k)
  q1b <- c(q1b, value)
}
plot(x= kappa, y = q1b, type = "l", lwd = 3, col = "dodgerblue4", 
     main = "MBS Price vs. Kappa", ylab = "PV0", xlab = "Kappa")
# c) for different r_bar
r_bar <- seq(0.03, 0.09, 0.01)
q1c <- c()
for (rbar in r_bar) {
  value <- findMBS(r_bar = rbar)
  q1c <- c(q1c, value)
}
plot(x= r_bar, y = q1c, type = "l", lwd = 3, col = "firebrick2", 
     main = "MBS Price vs. r_bar", ylab = "PV0", xlab = "r_bar")
# d) for different sigma
sig <- seq(0.1, 0.2, 0.01)
q1d <- c()
for (sigma in sig) {
  value <- findMBS(sig = sigma)
  q1d <- c(q1d, value)
}
plot(x= sig, y = q1d, type = "l", lwd = 3, col = "darkorchid4", 
     main = "MBS Price vs. Sigma", ylab = "PV0", xlab = "Sigma")


##### Question 2 #####
rPaths <- cirPath(r0=.078, sig=.12, kappa=.6, r_bar=.08,t = 30, nSim = 50000)
x <- uniroot(fitOAS, r_t = rPaths, lower = -0.015, upper = 0)$root

##### Question 3 #####
y <- 0.001
p0 <- 110000
p_plus <- fitOAS(x+y, rPaths, expPv = 0)
p_minus <- fitOAS(x-y, rPaths, expPv = 0)

duration <- (p_minus-p_plus)/(2*y*p0)
convexity <- (p_plus+p_minus-2*p0)/(2*p0*y^2)
