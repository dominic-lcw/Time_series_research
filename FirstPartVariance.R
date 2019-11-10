###-----------------------------------------------
### The data. **AR included**
###-----------------------------------------------
f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *10 + (t-0.2)^3 *20
}
ar = function(n, ar){ 
  z = rep(NA,n)
  z[1] = rnorm(1)
  for(i in 2:n){
    z[i] = ar*z[i-1] + rnorm(1,1)
  }
  z
}
n = 1000
ts.plot(ar(1000,0.8))

mu = f((1:n)/n)
z = c(rnorm(n/2,1), rnorm(n/2, 1,3)) 
x = z+mu
par(mfrow = c(1, 1))
ts.plot(x)

m1 = c(0.7071, -0.7071)
d1 = diff_seq(x, m1)^2

ts.plot(d1)
# abline(h = mean(d1[1:(n/2)]), col = 'red', lwd = 1.5)
# abline(h = mean(d1[(n/2+1):(n-1)]), col = 'blue', lwd = 1.5)
abline(h = mean(d1[1:n/2]), col = 'red', lwd = 1.5)
abline(h = mean(d1[(n/2+1):(n-1)]), col = 'blue', lwd = 1.5)



ts.plot(d1)
abline(h = mean(d1[1:k]), col = 'red', lwd = 1.5)
abline(h = mean(d1[k:(n-1)]), col = 'blue', lwd = 1.5)

b = KS_new(d1);b

k = K(d1, 0);k
k = 51
n = length(d1);n
d1_bar = mean(d1);x_bar

d2 = d1[1:k];
n1 = length(d2);n1
sigma_hat = var_head(d2, 2*n1^(1/3));sigma_hat

Tn =cumsum(d1 - d1_bar) / sqrt(n);ts.plot(Tn)
TS = max(abs(Tn / sqrt(sigma_hat)));TS


###----------------------------------
### Self Method
###----------------------------------
M <- function(x, t){
	#Standard KS
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

diff_seq = function(ts, diff){
	d = rep(NULL, length(ts))
	for (i in length(diff):length(ts)){
		d[i] = as.numeric(ts[i:(i-length(diff)+1)]%*%diff)
	}
	return(na.omit(d))
}

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

KS_new = function(x){ #
  k = K(x, 0)
  n = length(x)
  x_bar = mean(x)
  
  x1 = x[1:k]
  n1 = length(x1)
  sigma_hat = var_head(x1, 2*n1^(1/3))

  Tn =cumsum(x - x_bar) / sqrt(n)
  TS = max(abs(Tn / sqrt(sigma_hat)))
}

###-----------------------------------------------
### Normal Case. (Mid Point in the middle)
###-----------------------------------------------
m1 = c(0.7071, -0.7071)

if(1){
  n = 2000
  n_sim = 200
  delta = seq(from = 0, to = 1, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 2),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('self', 'self-new')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = 0
      #z = c(r(n/2,6), rt(n/2, 6)*(1+del))#random
      z = c(rnorm(n/2,1),rnorm(n/2,1)*(1+del))
      x = z+mu
      d1 = diff_seq(x, m1)^2
      out[i_sim, i_delta, 1] = KS_std(log(d1))
      out[i_sim, i_delta, 2] = KS_new(d1)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > 1.358
out0[,,2] = out[,,2] > 1.358
power = apply(out0, c(2,3), mean, na.rm = TRUE)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power n = 2000", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Self', 'Self new'), 
  col = c('blue','red'), lty = c(1,2,3), cex = 0.8)

#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)






