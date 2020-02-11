###-----------------------------------------------
### The data. **AR included**
###-----------------------------------------------
f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *10 + (t-0.2)^3 *20
}
ar = function(n, th1, th2){ 
  ar1 = rep(NA,n)
  ar1[1] = rnorm(1)
  ar1[2] = rnorm(1)
  ar1[3] = rnorm(1)
  for(i in 4:n){
    ar1[i] = th1*ar1[i-1] -th2*ar1[i-3]  + rnorm(1,1)
  }
  ar1
}
arma11 = function(n, th1, ph1){
	wn = rnorm(n)
	arma11  = rep(NA, n)
	arma11[1] = wn[1]
	for (i in 2:n){
		arma11[i] = arma11[i-1]*th1 + wn[i-1]*ph1 + wn[i]
	}
	arma11
}
ma = function(n, ph1, p){
  wn = rnorm(n)
  ma  = rep(NA, n)
  ma[1:p] = wn[1:p]
  for (i in (p+1):n){
    ma[i] = wn[(i-1):(i-p)]%*%ph1 + wn[i]
  }
  ma
}
n = 1600
#a = ar(1000,0.8, 0.3)
a = arma11(n, -0.8, -0.4)
mu = f((1:n)/n)
z = c(rnorm(n/2,1), rnorm(n/2, 1)*5) 
x = z+mu+a
par(mfrow = c(1, 1))
ts.plot(x)

###----------------------------------
### Self Method
###----------------------------------
diff_seq = function(ts, diff, h){
	d = rep(NULL, length(ts))
	m = length(diff)
	for (i in (m*h):length(ts)){
		d[i] = 0
		for (j in 1:m){
			d[i] = d[i] + ts[i-h*(j-1)]*diff[j]
		}
	}
	return(na.omit(d))
}

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

###----------------------------------
### Differencing Sequence
###----------------------------------
m1 = c(0.7071, -0.7071);sum(m1)
m2 = c(0.8090, -0.5, -0.3090);sum(m2)
m3 = c(0.1942, 0.2809, 0.3832, -0.8582);sum(m3)
m4 = c(0.2708, -0.0142, 0.6909, -0.4858, -0.4617);sum(m4)
m5 = c(0.9064, -0.2600, -0.2167, -0.1774, -0.1420, -0.1103);sum(m5)
m6 = c(0.2400, 0.0300, -0.0342, 0.7738, -0.3587, -0.3038, -0.3472);sum(m6)
m7 = c(0.9302, -0.1965, -0.1728, -0.1506, -0.1299, -0.1107, -0.0930, -0.0768);sum(m7)
m8 = c(0.2171, 0.0467, -0.0046, -0.0348, 0.8207, -0.2960, -0.2453, -0.2260, -0.2879);sum(m8)
m9 = c(0.9943, -0.1578, -0.1429, -0.1287, -0.1152, -0.1025, -0.0905, -0.0792, -0.0687,
 -0.0588);sum(m9)
m10 = c(0.1995, 0.0539, 0.0104, -0.0140, -0.0325, 0.8510, -0.2384, -0.2079, -0.1882,
 -0.1830,-0.2507);sum(m10)

###----------------------------------
### Simulation of different differencing sequences
###----------------------------------
mean_func <- function(n){
	x = f(1:n/n)
	return(x)
}

arma_func <- function(n){
	# a = ar(n, 0.7, 0.3)
	# a = arma11(n, -0.5, 0.5)
	a = ma(n, c(0.4,0.5,-0.6,-0.9,0,0,0,0.9), 8)
	return(a)
}

res <- function(n, del){
	z = c(rnorm(n/2, 0,1), rnorm(n/2, 0, 1+del))
	#z = c(rexp(n/2, 1)*1, rexp(n/2, 1)*(1+del)) #Exponential Distribution
	return(z)
}

if(1){
	n = 400
	n_sim = 2^9
	order =5
	delta = seq(from = 0, to = 2, length.out = 21)
	out = array(NA, dim = c(n_sim, length(delta), order),
				dimnames = list(paste0('isim=', 1:n_sim),
								paste0('delta=', delta),
								paste0('m=',1:order)))
	for(i_sim in 1:n_sim){
		set.seed(i_sim)
		for(i_delta in 1:length(delta)){
			del = delta[i_delta]
			mu = mean_func(n)
			#a = arma_func(n)
			z = res(n, del)
			#x = z+mu+a
			x = z+mu
			d1 = diff_seq(x, m1,1)
			d2 = diff_seq(x, m2,1)
			d3 = diff_seq(x, m3,1)
			d5 = diff_seq(x, m5,1)
			d9 = diff_seq(x, m9,1)
			out[i_sim, i_delta, 1] = KS_std(log(d1^2))
			out[i_sim, i_delta, 2] = KS_std(log(d2^2))
			out[i_sim, i_delta, 3] = KS_std(log(d3^2))
			out[i_sim, i_delta, 4] = KS_std(log(d5^2))
			out[i_sim, i_delta, 5] = KS_std(log(d9^2))
			
		}
		if(i_sim%%20==0) cat(i_sim, ' >> ')
	}
}

###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > 1.358
out0[,,2] = out[,,2] > 1.358
out0[,,3] = out[,,3] > 1.358
out0[,,4] = out[,,4] > 1.358
out0[,,5] = out[,,5] > 1.358
power = apply(out0, c(2,3), mean)
power
###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 1))
#Plot 1
matplot(delta, 100*power, type = 'b', col = 1:10, lty = c(1,2,3), 
	pch = "12359", main = "Power (added ARMA series)", ylim = c(0, 100),
    ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
# legend('bottomright', legend = c(paste0("m=",c(1,)), 
#   col = 1:10, lty = c(1,2,3), cex = 0.4)

###-----------------------------------------------------------
### Free Space
###-----------------------------------------------------------
### +ve acf lag1
n = 1000
a = arma11(n, 0.8, -0.4)
z = c(rnorm(n/2,1), rnorm(n/2, 1)*3) 
x = z+a
par(mfrow = c(1,1))
ts.plot(x)
acf(x)

m1 = c(1, -1)
d1 = diff_seq(x, m1)
m2 = c(1, -2, 1)
d2 = diff_seq(x, m2)

par(mfrow = c(2, 1))
ts.plot(d1)
ts.plot(d2)
var(d1)
var(d2)

### -ve acf lag1
n = 1000
a = arma11(n, -0.8, -0.4)
z = c(rnorm(n/2,1), rnorm(n/2, 1)*3) 
x = z+a
par(mfrow = c(1,1))
ts.plot(x)
acf(x)

m1 = c(1, -1)
d1 = diff_seq(x, m1)
m2 = c(1, -2, 1)
d2 = diff_seq(x, m2)

par(mfrow = c(2, 1))
ts.plot(d1)
ts.plot(d2)
var(d1)
var(d2)

###Compare 
sum(rexp(9, 1));rgamma(1, 9, 1)
set.seed(111)
hist(m1[1]*rnorm(1000) - m1[2]*rnorm(1000))
hist(m2[1]*rnorm(1000) -m2[2]*rnorm(1000) + m2[3]*rnorm(1000))


m1 = c(0.7071, -0.7071);sum(m1)
m2 = c(0.8090, -0.5, -0.3090);sum(m2)


sum(m1^2)
sum(m2^2)


###Diff 18
n = length(x)
x_bar = mean(x)
sigma_hat = var_head(x, 2*n^(1/3))
Tn =cumsum(x - x_bar) / sqrt(n)
TS = max(abs(Tn / sqrt(sigma_hat)))
ts.plot(d18)

x = log(d18^2)
ts.plot(Tn)

TS

ts.plot(x)




