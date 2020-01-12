###-----------------------------------------------
### The data. **AR included**
###-----------------------------------------------
f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *10 + (t-0.2)^3 *20
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
z = c(rnorm(n/2,3), rnorm(n/2, 3)*3) 
x = z+mu
par(mfrow = c(1, 1))
ts.plot(x)
d = diff(x)/sqrt(2)
d2 = d^2
ts.plot(d)
ts.plot(d2)
abline(h = mean(d2[1:(n/2)]), col = 'red', lwd = 1.5)
abline(h = mean(d2[(n/2+1):(n-1)]), col = 'blue', lwd = 1.5)


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

###-----------------------------------------------
### Textbook method
###-----------------------------------------------
M <- function(x, t){
  n = length(x)
  bar_x = mean(x)
  (sum((x[1:t] - bar_x)^2) - t/n*sum((x - bar_x)^2)) / (sqrt(n))
}

K <- function(x, gam){
  n = length(x)
  k_t = rep(NULL, (n-1))
  for (i in 1:(n-1)){
    k_t[i] = M(x, i)/((i/(n-1))*(1-i/(n-1)))^gam
  }
  which.min(k_t)
}

V_k <- function(k, x){
  bar = mean(x)
  n = length(x)
  (sum((x-bar)^4) - (1/k)*(sum((x[1:k]-bar)^2)^2) - (1/(n-k))*(sum((x[(k+1):n]-bar)^2, na.rm = TRUE)^2))/n
}

textbook <- function(d){
  k = K(d, 0)
  v = V_k(k, d)
  M_t <- rep(NULL, n)
  for (i in 1:(n-1)){
    M_t[i] = M(d, i) / sqrt(v)
  }
  max(abs(M_t))
}

###-----------------------------------------------
### Liver Method
###-----------------------------------------------
L = function(y, g, tau){
  n = length(y)
  l1 = tau * log(1/tau * sum((y[1:tau] - g[1:tau])^2))
  if (tau !=n){
    l2 = (n - tau) * log(1/(n-tau)*sum((y[(tau+1):n] - g[(tau+1):n])^2))}
  else{l2 = 0}
  l = l1 + l2
}

MaxL = function(y, g){
  n = length(y)
  max_d= 0
  max_i = 0
  Ln = L(y, g, n)
  for (i in 1:n){
    d = Ln -  L(y, g, i)
    if(d > max_d){max_d = d; max_i = i}
  }
  c(max_d, max_i)
}

MLsigma1 = function(y, g, tau){sigma = 1/tau * sum((y[1:tau] - g[1:tau])^2)}
MLsigma2 = function(y, g, tau){n = length(y);sigma = 1/(n-tau)*sum((y[(tau+1):n] - g[(tau+1):n])^2)}

Liver_test = function(d, alpha){
  a = 2*log(log(n))^(1/2)/log(n)
  b = (2*log(log(n)) + 1/2*log(log(log(n))) - log(gamma(1/2)))/log(n)
  T = a*sqrt(log(n))*sqrt(d) - b*log(n)
}

spline_test = function(x){
  epsilon = 0.00001
  n = length(x)
  ax = (1:n)/n
  old_d = 0
  count = 0
  weight = rep(1:n)
  while (count<50){
    g = smooth.spline(ax, x, w = weight, cv = TRUE)$y
    result= MaxL(x, g)
    new_d = result[1]
    tau = result[2]
    if ((new_d - old_d)<epsilon){
      return(Liver_test(new_d, 0.05))
    }else{
      sigma1 = MLsigma1(x, g, tau)
      sigma2 = MLsigma2(x, g, tau)
      weight = c(rep((1/sigma1)^2, n/2), rep((1/sigma2)^2, n/2))
    }
    count = count+1
    old_d = new_d  
  }
  return(Liver_test(new_d, 0.05))
}
###-----------------------------------------------
### Normal Case. (Mid Point in the middle)
###-----------------------------------------------

del = 0.1
if(1){
  n = 400
  n_sim = 200
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('liver', 'textbook','log(d^2)','d^2')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = f((1:n)/n)
      z = c(rnorm(n/4, 0, 1), rnorm(n*3/4, 0, 1)*(1+del))
      #z = c(r(n/2,6), rt(n/2, 6)*(1+del))#random
      # z = c(rnorm(n/4,1),rnorm(n/4,1)*(2+del), rnorm(n/2,1)*(7-del))
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = spline_test(x)
      out[i_sim, i_delta, 2] = textbook(x)
      out[i_sim, i_delta, 3] = KS_std(log(d^2))
      out[i_sim, i_delta, 4] = KS_std(d^2)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out[,,2] > 1.358
out0[,,3] = out[,,3] > 1.358
out0[,,4] = out[,,4] > 1.358
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (n_sim=1000, n=500)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Liver', 'Textbook','log(d^2)','d^2'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.4)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)

###-----------------------------------------------
### Changing Mean
###-----------------------------------------------

if(1){
  n = 400
  n_sim = 200
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('liver', 'textbook','log(d^2)','d^2')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = f((1:n)/n)
      #z = c(r(n/2,6), rt(n/2, 6)*(1+del))#random
      z = c(rnorm(n/2,1),rnorm(n/2,1)*(1+del))
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = spline_test(x)
      out[i_sim, i_delta, 2] = textbook(x)
      out[i_sim, i_delta, 3] = KS_std(log(d^2))
      out[i_sim, i_delta, 4] = KS_std(d^2)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out[,,2] > 1.358
out0[,,3] = out[,,3] > 1.358
out0[,,4] = out[,,4] > 1.358
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (n_sim=1000, n=500)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Liver', 'Textbook','log(d^2)','d^2'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.4)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)

###-----------------------------------------------
### Early Change Point
###-----------------------------------------------

if(1){
  n = 400
  n_sim = 200
  delta = seq(from = 0, to = 5, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('liver', 'textbook','log(d^2)','d^2')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = 0
      #z = c(r(n/2,6), rt(n/2, 6)*(1+del))#random
      #z = c(rnorm(n*2/8,1),rnorm(n*7/8,1)*(1+del))
      # z = c(rexp(n*1/2, 1), rexp(n*1/2, 1)*(1+del))
      p1 = c(rbinom(n*1/8, 1, 0.5));p2 = c(rbinom(n*7/8, 1, 0.5));z = c(rnorm(n*1/8,-2)*p1+rnorm(n*1/8,2)*(1-p1),(rnorm(n*7/8,-2)*p2+rnorm(n*7/8,2)*(1-p2))*(1+del))
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = spline_test(x)
      out[i_sim, i_delta, 2] = textbook(x)
      out[i_sim, i_delta, 3] = KS_std(log(d^2))
      out[i_sim, i_delta, 4] = KS_std(d^2)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

#Test of multinomial case
del = 0.5
n = 400
p1 = c(rbinom(n*1/2, 1, 0.5))
p2 = c(rbinom(n*1/2, 1, 0.5))
z = c(rnorm(n*1/2,-2*p1+rnorm(n*1/2,10)*(1-p1),(rnorm(n*1/2,-10)*p2+rnorm(n*1/2,10)*(1-p2))*(1+del))
ts.plot(z)



###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out[,,2] > 1.358
out0[,,3] = out[,,3] > 1.358
out0[,,4] = out[,,4] > 1.358
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Liver', 'Textbook','log(d^2)','d^2'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)

###-----------------------------------------------
### Early Change Point
###-----------------------------------------------

if(1){
  n = 400
  n_sim = 200
  delta = seq(from = 0, to = 1, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 3),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('liver', 'textbook','log(d^2)')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = 0
      #z = c(r(n/2,6), rt(n/2, 6)*(1+del))#random
      #z = c(rnorm(n*2/8,1),rnorm(n*7/8,1)*(1+del))
      #z = c(rexp(n*1/2, 1), rexp(n*1/2, 1)*(1+del))
      z_bi = c(rnorm(n*1/2,-2,1)*p1+rnorm(n*1/2,10,1)*(1-p1),
      (rnorm(n*1/2,-10,(1+del))*p2+rnorm(n*1/2,10,(1+del))*(1-p2))) #Bimodal T-distribution
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = spline_test(x)
      out[i_sim, i_delta, 2] = textbook(x)
      out[i_sim, i_delta, 3] = KS_std(log(d^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out[,,2] > 1.358
out0[,,3] = out[,,3] > 1.358
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen')
#Plot 1
matplot(delta, 100*power, type = 'b', col = col, pch= 1, lty = c(1,2,3),
     main = "Power", ylim = c(0, 100), ylab = expression(K(delta)~"/%"), 
     xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Liver', 'Textbook','log(d^2)'), 
  col = c('blue','red','darkgreen'), lty = c(1,2,3), cex = 0.6)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)

