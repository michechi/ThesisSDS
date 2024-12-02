# IMPLEMENTATION JULIA CODE AVAILABLE AT https://github.com/jwmi/BayesianMixtures.jl

logsumexp <- function(a,b){
  m = max(a,b)
  if(m == -Inf){return(-Inf)}else{return(log(exp(a-m) + exp(b-m)) + m)}
}

# Compute log_v[t] = log(V_n(t)) under the given MFM parameters, for t=1:upto.

log_Vn.M <- function(gamma_, n, upto, par){
  tolerance <- 1e-100
  log_v <- rep(0, upto)
  for(t in 1:upto){
    if(t >n){log_v[t] = -Inf}else{
      a = 0
      c = -Inf
      k = 1
      p = 0 
      dif = a-c
      #print(p < (1 - tolerance))
      #print(abs(a-c)> tolerance)
      
      while((abs(dif) > tolerance) || (p < (1 - tolerance))){ # Note: The first condition is false when a = c = -Inf
        if(k >= t){
          a = c
          b = lgamma(k+1)-lgamma(k-t+1)-lgamma(k*gamma_+n)+lgamma(k*gamma_)+dtpois(k,par,0, log = T)
          c = logsumexp(a,b)
        }
        p = p+exp(log(dtpois(k,1,0)))
        k = k+1
        if(a == c){dif = 0}else{dif = a-c}
        #print(p < (1 - tolerance))
        #print(abs(a-c)> tolerance)
      }
      log_v[t] = c
    }
  }
  return(log_v)
}


