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
ar(1000,0.8)

mu = f((1:n)/n)
z = c(rnorm(n/2,1), rnorm(n/2, 1)*3) 
x = z+mu
par(mfrow = c(1, 1))
ts.plot(x)

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
  k = K(x, 0)
  n = length(x)
  x_bar = mean(x)
  
  x1 = x[1:k]
  n1 = length(x1)
  sigma_hat = var_head(x1, 2*n1^(1/3))

  Tn =cumsum(x - x_bar) / sqrt(n)
  TS = max(abs(Tn / sqrt(sigma_hat)))
}







