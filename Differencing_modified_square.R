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


d = diff(x)/sqrt(2)
d2 = d^2
ts.plot(d)
ts.plot(d2)
abline(h = mean(d2[1:500]), col = 'red', lwd = 1.5)
abline(h = mean(d2[501:999]), col = 'blue', lwd = 1.5)
###----------------------------------
### Self Method
###----------------------------------

k_split = function(x){
  n = length(x)
  x_bar = mean(x)
  k = which.min(cumsum(x - x_bar)[1:(n-1)])
}

two_sigma = function(x,k){ #CP in mean and variance together
  n = length(x)
  x_bar1 = mean(x[1:k])
  x_bar2 = mean(x[(k+1):n]) 
  sigma1 = 1/k * sum((x[1:k] - x_bar1)^2)
  sigma2 = 1/(n-k) * sum((x[(k+1):n] - x_bar2)^2)
  c(sigma1, sigma2)
}

KS_std_new = function(d){ #
  k = k_split(d^2)
  sd = two_sigma(d,k)
  n = length(d)
  sd_seq = sqrt(c(rep(sd[1],k),rep(sd[2],(n-k))))
  dd = d/sd_seq + sd_seq
  sigma_hat = var_head(dd, 2*n^(1/3))
  dd_bar = mean(dd)
  Tn =cumsum(dd - dd_bar) / sqrt(n)
  max(abs(Tn / sqrt(sigma_hat)))
}

k = k_split(d2);k
two_sd = two_sigma(d, k);two_sd
sd_seq = sqrt(c(rep(two_sd[1],k),rep(two_sd[2],(n-k))));sd_seq
dd = d/sd_seq + sd_seq;ts.plot(dd)
n = length(dd);n
sigma_hat = var_head(dd, 2*n^(1/3));sigma_hat
dd_bar = mean(dd)
Tn =cumsum(dd - dd_bar) / sqrt(n);ts.plot(Tn)
max(abs(Tn / sqrt(sigma_hat)))





###-------------------------------------
### Simulation
###-------------------------------------
if(1){
  n = 400
  n_sim = 200
  delta = seq(from = 0, to = 1, by = 0.1)
  out = array(NA, dim =c(n_sim, length(delta), 1),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('Transfer')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = 0
      #z = c(rt(n/2,6), rt(n/2, 6)*(1+del))#random
      z = c(rnorm(n/2,1),rnorm(n/2,1)*(1+del))
      x = z+mu
      d = diff(x)/sqrt(2)
      out[i_sim, i_delta, 1] = KS_std_new(d^2)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > 1.358
power = apply(out0, c(2,3), mean, na.rm = TRUE)
power


###-------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'b', col = col, pch = 1, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('New Method'), 
  col = c('blue'), lty = c(1,2,3), cex = 0.8)