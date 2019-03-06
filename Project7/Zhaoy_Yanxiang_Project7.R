# Zhao_Yanxiang_Projcet7
source("finiteDiff.r")
source("bsm.r")
########################### Main ########################### 
##### Question 1 #####
# set parameters
s0 <- 10
k <- 10
r <- 0.04
sig <- 0.2
t <- 0.5
dt <- 0.002
dx <- sig*sqrt(dt)

# i. Values
# (a) Explicit Finite-Difference Method
Pa <- efd(type="put", euro=T, s0, k, r, sig, t, dt, dx, log=T)
# (b) Implicit Finite-Difference Method
Pb <- ifd(type="put", euro=T, s0, k, r, sig, t, dt, dx, log=T)
# (c) Crank-Nicolson Finite-Difference Method
Pc <- cnfd(type="put", euro=T, s0, k, r, sig, t, dt, dx, log=T)
# output
c(EFD = Pa, IFD = Pb, CNFE = Pc)

# ii. Comparison
s0 <- 4:16
p_edf <- p_idf <- p_cndf <- p_bs <- vector()
for (i in 1:length(s0)){
  p_edf[i] <- efd(type="put", euro=T, s0[i], k, r, sig, t, dt, dx, log=T)
  p_idf[i] <- ifd(type="put", euro=T, s0[i], k, r, sig, t, dt, dx, log=T)
  p_cndf[i] <- cnfd(type="put", euro=T, s0[i], k, r, sig, t, dt, dx, log=T)
  p_bs[i] <- bsPut(s0[i], t, k, r, sig)
}
err <- cbind(p_edf, p_idf, p_cndf)-p_bs
rownames(err) <- s0
# plot
matplot(s0, abs(err), type ="l", lwd = 3, lty=1, main = "Absolute error against BSM Model")
legend("topleft", legend = colnames(err), lwd = 3, col = 1:3)

##### Question 2 #####
# set parameters
s0 <- 10
k <- 10
r <- 0.04
sig <- 0.2
t <- 0.5
dt <- 0.002
ds <- 0.25

# i. Values
# (a) Explicit Finite-Difference Method
Ca <- efd(type="call", euro=F, s0, k, r, sig, t, dt, ds, log=F)
Pa <- efd(type="put", euro=F, s0, k, r, sig, t, dt, ds, log=F)
# (b) Implicit Finite-Difference Method
Cb <- ifd(type="call", euro=F, s0, k, r, sig, t, dt, ds, log=F)
Pb <- ifd(type="put", euro=F, s0, k, r, sig, t, dt, ds, log=F)
# (c) Crank-Nicolson Finite-Difference Method
Cc <- cnfd(type="call", euro=F, s0, k, r, sig, t, dt, ds, log=F)
Pc <- cnfd(type="put", euro=F, s0, k, r, sig, t, dt, ds, log=F)
# output
cbind(EFD = c(call = Ca, put = Pa), IFD = c(Cb, Pb), CNFE = c(Cc, Pc))

# ii. Graphs
s0 <- 4:16
c_edf2 <- c_idf2 <- c_cndf2 <- vector()
p_edf2 <- p_idf2 <- p_cndf2 <- vector()
for (i in 1:length(s0)){
  # calls
  c_edf2[i] <- efd(type="call", euro=F, s0[i], k, r, sig, t, dt, ds, log=F)
  c_idf2[i] <- ifd(type="call", euro=F, s0[i], k, r, sig, t, dt, ds, log=F)
  c_cndf2[i] <- cnfd(type="call", euro=F, s0[i], k, r, sig, t, dt, ds, log=F)
  # puts
  p_edf2[i] <- efd(type="put", euro=F, s0[i], k, r, sig, t, dt, ds, log=F)
  p_idf2[i] <- ifd(type="put", euro=F, s0[i], k, r, sig, t, dt, ds, log=F)
  p_cndf2[i] <- cnfd(type="put", euro=F, s0[i], k, r, sig, t, dt, ds, log=F)
}
matplot(s0, cbind(c_edf2, c_idf2, c_cndf2), type ="l", lwd = 1, lty=1,
        main = "Call Prices", xlab = "Stock price", ylab = "Call Price")
legend("topleft", legend = c("EDF", "IDF", "CNDF"), lwd = 1, col = 1:3)
matplot(s0, cbind(p_edf2, p_idf2, p_cndf2), type ="l", lwd = 1, lty=1,
        main = "Put Prices", xlab = "Stock price", ylab = "Put Price")
legend("topright", legend = c("EDF", "IDF", "CNDF"), lwd = 1, col = 1:3)
