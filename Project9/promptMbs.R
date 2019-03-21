source('abs.R')

print("please enter the following inputs for the MBS: ")
pv0 <- as.numeric(readline("The notional amount of the loan is ($): ")) 
wac <- as.numeric(readline("The weighted average coupon is (%): "))/100
years <- as.numeric(readline("The maturity of the loan is (years): ")) 
r0 <- as.numeric(readline("The initial interest rate is (%): "))/100
r_bar <- as.numeric(readline("The average interest rate is (%): "))/100
sig <- as.numeric(readline("The annual volatility is (%): "))/100
kappa <- as.numeric(readline("The mean reversion coefficient is: "))
nSim <- as.numeric(readline("The number of simulations to run: "))
cat('The price of the MBS is: $')
cat(findMBS(pv0, r0, r_bar, wac, kappa, sig, years, nSim))