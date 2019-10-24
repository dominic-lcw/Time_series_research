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
z = c(rt(n/2,3), rt(n/2, 3)*4) 
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
  k_t = rep(NULL, n)
  for (i in 1:n){
    k_t[i] = M(x, i)/((i/n)*(1-i/n))^gam
  }
  which.min(k_t)
}

V_k <- function(k, x){
  bar = mean(x)
  n = length(x)
  (sum((x-bar)^4) - (1/k)*(sum((x[1:k]-bar)^2)^2) - (1/(n-k))*(sum((x[(k+1):n]-bar)^2, na.rm = TRUE)^2))/n
}

textbook <- function(d){
  k = K(d, 1/3)
  v = V_k(k, d)
  M_t <- rep(NULL, n)
  for (i in 1:(n-1)){
    M_t[i] = M(d, i) / sqrt(v)
  }
  max(abs(M_t))
}

#Testing of code
k = K(d, 1/3)
v = V_k(k, d)
M_t <- rep(NULL, n)
for (i in 1:(n-1)){
  M_t[i] = M(d, i) / sqrt(v)
}
ts.plot(M_t)

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
### Simulation Part
###-----------------------------------------------

if(1){
  n = 500
  n_sim = 1000
  delta = seq(from = 0, to = 1, by = 0.05)
  out = array(NA, dim =c(n_sim, length(delta), 3),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('liver', 'textbook','self')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = f((1:n/n))
      z = c(rt(n/2,6), rt(n/2, 6)*(1+del))
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = spline_test(x)
      out[i_sim, i_delta, 2] = textbook(d)
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
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend(0.75, 20, legend = c('Liver', 'Textbook','self'), col = c('blue','red','darkgreen'), lty = c(1,2,3), cex = 0.5)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)

