###-----------------------------------------------------------
### Change in mean (KS test)
###-----------------------------------------------------------
set.seed(1)
f = function(){
  rep(c(1, 2), c(n/2, n/2))
}
n = 1000
mu = f()
z = c(rnorm(n))
x = z+mu
ts.plot(x)

###-----------------------------------------------------------
### Standard KS test-statistics
###-----------------------------------------------------------
kern = function(x){ #Bartlett kernel
  ifelse(abs(x)<=1,1-abs(x), 0)
}

var_head = function(x, l){ #Kernel Estimation of volatility
  n = length(x)
  l = floor(l)
  sum = 0
  covar = c(acf(x, l+1, type = 'cov', plot=FALSE)$acf)
  k = (-l):l
  out = sum(kern(k/l) * covar[abs(k) +1])
  
  #for (k in -l:l)
  #  sum = sum + kern(k/l) * covar[abs(k) +1] 
  #}
  return(out)
}

KS_std = function(x){ #
  n = length(x)
  x_bar = mean(x)
  sigma_hat = var_head(x, 2*n^(1/3))
  Tn =cumsum(x - x_bar) / sqrt(n)
  TS = max(abs(Tn / sqrt(sigma_hat)))
}

###-----------------------------------------------------------
### SN KS test-statistics
###-----------------------------------------------------------
KS_sn = function(x){
  n = length(x)
  x_bar = mean(x)
  Tn =cumsum(x - x_bar) / sqrt(n)
  V = rep(NA, n-1)
  
  s = cumsum(x)
  rs =  rev(cumsum(rev(x))) #sum backward
  
  for (k in 1:(n-1)){
    t = 1:k
    rt = (k+1):n
    V[k] = sum((s[t] - (t/k) * s[k]) ^2) + sum((rs[rt] - (n-rt+1) / (n-k) * rs[k+1])^2)
    #V[k] = (sum((s[t] - (t/k)*s[k])^2) + sum((rs[rt] - (n-rt+1)/(n-k) * rs[k+1])^2))
  }
  V = V/n^2
  TS = max(abs(Tn[-n]^2 / V))
}

###-----------------------------------------------------------
### Simulation Experiments
###-----------------------------------------------------------

if(1){
  n = 1000
  n_sim = 100
  delta = seq(from = 0, to = 2, by = .05)
  out = array(NA, dim = c(n_sim, length(delta), 2),
              dimnames = list(paste0('iSim=', 1:n_sim), paste0("delta=", delta), c("std", 'sn')))
  for(i_sim in 1:n_sim){
    for(i_delta in 1:length(delta)){
    set.seed(i_sim)
      d = delta[i_delta]
      x = rnorm(n) + rep(c(0, d), c(n/2, n/2))
      out[i_sim, i_delta, 1] = KS_std(x)
      out[i_sim, i_delta, 2] = KS_sn(x)
    }
    if(i_sim%%10 ==0) cat(i_sim, " >> ")
  }
}
###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > 1.358
out0[,,2] = out[,,2] > 40.1
power = apply(out0, c(2,3), mean)

###-----------------------------------------------------------
### Results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue', 'red')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.2) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.15))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)



