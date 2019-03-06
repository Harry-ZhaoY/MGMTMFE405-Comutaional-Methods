# Zhao_Yanxiang_Projcet8
setwd("C:/Users/harry/OneDrive/Documents/GitHub/MGMTMFE405-Comutaional-Methods/Project8")
source('fixedIncome.R')
########################### MAIN ########################### 
##### Question 1 #####
# set parameters
r0 <- .05
sig <- .18
kappa <- .82
r_bar <- .05
nSim <- 10000
fv <- 1000
# (a) pure discout bond
t <- 0.5
# generate paths
paths1a <- abs(vasicekPath(r0, sig, kappa, r_bar, t, nSim))
# find pure discount bond price
bond_1a <- zcBond(paths1a, fv)
# (b) coupon paying bond
t <- 4
c <- 30
# generate paths
paths1b <- abs(vasicekPath(r0, sig, kappa, r_bar, t, nSim))
# find pure discount bond price
bond_1b <- cpnBond(paths1b, t, c, fv)
# (c) European Call option on pure discount bond
exerT <- .25
t <- .5
k <- 980
fv <- 1000
zcOption_vasicek(paths1a, sig, kappa, r_bar, exerT, t, k, fv, 'call')
# (d) European Call option on coupon bond

