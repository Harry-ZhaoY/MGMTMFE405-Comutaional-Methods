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
nSim <- 1000
fv <- 1000
# (a) pure discout bond
t <- 0.5
# generate paths
paths1a <- vasicekPath(r0, sig, kappa, r_bar, t, nSim)
# find pure discount bond price
bond_1a <- zcBond(paths1a, fv)
# (b) coupon paying bond
t <- 4
pmtT <- seq(0.5, t, 0.5)
c <- 30
# generate paths
paths1b <- vasicekPath(r0, sig, kappa, r_bar, t, nSim)
# find pure discount bond price
bond_1b <- cpnBond(paths1b, pmtT, c, fv)
# (c) European Call option on pure discount bond
exerT <- .25
t <- .5
k <- 980
opt_1c <- zcOption_vasicek(paths1a, sig, kappa, r_bar, exerT, t, k, fv, 'call')
# (d) European Call option on coupon bond
exerT <- .25
t <- 4
k <- 980
c <- 30
mSim <- 100
paths1d <- vasicekPath(r0, sig, kappa, r_bar, exerT, nSim)
opt_1d <- cpnOption(paths1d, sig, kappa, r_bar, exerT, t, k, c, fv, "call", "Vasicek", mSim)

##### Question 2 #####
# set parameters
r0 <- .05
sig <- .18
kappa <- .92
r_bar <- .055
nSim <- 1000
fv <- 1000
# (a) European Call option
exerT <- 0.5
t <- 1
k <- 980
mSim <- 100
paths2a <- cirPath(r0, sig, kappa, r_bar, t, nSim)
opt_2a <- zcOption(paths2a, sig, kappa, r_bar, exerT, t, k, fv, "call", "CIR", mSim)
# (b) Explicite solution
opt_2b <- zcCall_cir(r0, sig, kappa, r_bar, exerT, t, k, fv)
##### Question 3 #####
# set parameters
x0 <- y0 <- 0
r0 <- phi_t <- 0.03
rho <- 0.7
a <- 0.1
b <- 0.3
sig <- 0.03
eta <- 0.08
t <- 0.5
nSim <- 1000
# paths
paths3 <- g2ppPath(x0, y0, r0, phi_t, rho, a, b, sig, eta, t, nSim)
# price the put option
exerT <- 0.5
t <- 1
k <- 985
mSim <- 100
opt_3 <- zcOption_gcpp(paths3, x0, y0, r0, phi_t, rho, a, b, sig, eta, exerT, t, k, fv, "put", mSim)
