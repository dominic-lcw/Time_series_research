###-----------------------------------------------
### The data. **AR included**
###-----------------------------------------------
f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *20 + (t-0.2)^3 *40
}
ar = function(n, ar, df){
  z = rep(NA,n)
  z[1] = rt(1,df)
  for(i in 2:n){
    z[i] = ar*z[i-1] + rt(1,df)
  }
  z
}
n = 1000
mu = f((1:n)/n)
z = c(rnorm(n/2), rnorm(n/2))
#z = c(rt(n/2,3), rt(n/2, 3)) 
x = z+mu
par(mfrow = c(1, 1))
ts.plot(x)

###----------------------------------
### Self Method
###----------------------------------

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
  return(out)
}

KS_std = function(x){ #
  n = length(x)
  x_bar = mean(x)
  sigma_hat = var_head(x, 2*n^(1/3))
  Tn =cumsum(x - x_bar) / sqrt(n)
  TS = max(abs(Tn / sqrt(sigma_hat)))
}

###-------------------------------------
### Simulation
###-------------------------------------
if(1){
  n = 400
  n_sim = 1000
  delta = seq(from = 0, to = 1, by = 0.1)
  out = array(NA, dim =c(n_sim, length(delta), 2),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('X^2', 'Log(X^2)')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = 0
      z = c(rt(n/2,6), rt(n/2, 6)*(1+del))#random
      #z = c(rnorm(n/2,1),rnorm(n/2,1)*(1+del))
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = KS_std(d^2)
      out[i_sim, i_delta, 2] = KS_std(log(d^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

###-------------------------------------------------------
### Plot 
###-------------------------------------------------------
mu = 0
z = c(rt(n/2,6), rt(n/2, 6)*(1+del))#random
#z = c(rnorm(n/2,1),rnorm(n/2,1)*(1+del))
x = z+mu
d = diff(x)/sqrt(2)
par(mfrow = c(1,1))
ts.plot(d^2)
ts.plot(log(d^2))

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > 1.398
out0[,,2] = out[,,2] > 1.398
power = apply(out0, c(2,3), mean)
power
###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (n_sim=1000, n=500)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend(0.7, 25, legend = c('X^2', 'Log(X^2)'), col = c('blue','red'), 
       lty = c(1,2,3), cex = 0.6)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)



