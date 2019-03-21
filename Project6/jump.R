jump = function(lambda1=0.2, lambda2=0.4, capT=5){
  V0 = 20000
  L0 = 22000
  mu = -0.1
  sigma = 0.2
  gamma = -0.4
  r0 = 0.02
  n = capT*12
  delta = 0.25
  alpha = 0.7
  epsilon = 0.95
  beta = (epsilon-alpha)/capT
  r = (r0+delta*lambda2)/12
  pmt = L0*r/(1-1/(1+r)^n)
  a = pmt/r
  b = pmt/(r*(1+r)^n)
  c = 1+r
  steps3 = 100
  paths3 = 10000
  dt3 = capT/steps3
  V = vector("numeric")
  S = firstJump = 0
  V[1] = V0
  loan_Q = value_Q = value_S = 0
  default_option = default_count = expected_time = 0
  for(i in 1:paths3){
    payoff = 0
    tau = Q = capT
    # time of the first adverse
    S = which(rpois(steps3, lambda2*dt3) == 1)[1]
    if(is.na(S)){
      S = steps3+1
    }
    S= S*dt3
    for(j in 1:steps3){
      t = (j-1)*dt3
      Lt = a - b*c^(12*t)
      qt = alpha + beta*t
      # when adverse event happened
      if(S < t){
        value_S = V[j]
        tau = S
        break
      }
      # when collaertized value less than loan
      if(V[j] <= qt*Lt){
        Q = t
        loan_Q = Lt
        value_Q = V[j]
        tau = Q
        break
      }
      V[j+1] = V[j] + V[j]*(mu*dt3 + sigma*sqrt(dt3)*rnorm(1,0,1) +
                              gamma*rpois(1, lambda1*dt3))
    }
    # when tau is less than T, then exercise the option
    if(tau < capT){
      default_count = default_count + 1
      expected_time = expected_time + tau
      if(Q < S){
        payoff = max(0, loan_Q - epsilon*value_Q)
        default_option = default_option + payoff*exp(-r0*Q)
      }
      else{
        payoff = abs(a - b*c^(12*S) - epsilon*value_S)
        default_option = default_option + payoff*exp(-r0*S)
      }
    }
  }
  result = vector()
  result[1] = default_option/paths3
  result[2] = default_count/paths3
  result[3] = expected_time/default_count
  return(result)
}

      