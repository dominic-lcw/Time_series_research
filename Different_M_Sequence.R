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
z = c(rnorm(n/2,1), rnorm(n/2, 1)*3) 
x = z+mu
par(mfrow = c(1, 1))
ts.plot(x)


diff_seq = function(diff, ts){
	d = rep(NULL, length(ts))
	for (i in length(diff):length(ts)){
		d[i] = as.numeric(ts[i:(i-length(diff)+1)]%*%diff)
	}
	return(na.omit(d))
}

ds = c(1,-1)/sqrt(2)
d1 = diff_seq(ds, x)
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

if(1){
	n = 400
	n_sim = 5000
	order =10
	delta = seq(from = 0, to = 1, length.out = 21)
	out = array(NA, dim = c(n_sim, length(delta), order),
				dimnames = list(paste0('isim=', 1:n_sim),
								paste0('delta=', delta),
								paste0('m=',1:order)))
	for(i_sim in 1:n_sim){
		set.seed(i_sim)
		for(i_delta in 1:length(delta)){
			del = delta[i_delta]
			mu = 0
			z = c(rnorm(n/2, 1), rnorm(n/2, 1)*(1+del))
			x = z+mu
			d1 = diff_seq(m1, x)^2
			d2 = diff_seq(m2, x)^2
			d3 = diff_seq(m3, x)^2
			d4 = diff_seq(m4, x)^2
			d5 = diff_seq(m5, x)^2
			d6 = diff_seq(m6, x)^2
			d7 = diff_seq(m7, x)^2
			d8 = diff_seq(m8, x)^2
			d9 = diff_seq(m9, x)^2
			d10 = diff_seq(m10, x)^2
			out[i_sim, i_delta, 1] = KS_std(log(d1^2))
			out[i_sim, i_delta, 2] = KS_std(log(d2^2))
			out[i_sim, i_delta, 3] = KS_std(log(d3^2))
			out[i_sim, i_delta, 4] = KS_std(log(d4^2))
			out[i_sim, i_delta, 5] = KS_std(log(d5^2))
			out[i_sim, i_delta, 6] = KS_std(log(d6^2))
			out[i_sim, i_delta, 7] = KS_std(log(d7^2))
			out[i_sim, i_delta, 8] = KS_std(log(d8^2))
			out[i_sim, i_delta, 9] = KS_std(log(d9^2))
			out[i_sim, i_delta, 10] = KS_std(log(d10^2))
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
out0[,,6] = out[,,6] > 1.358
out0[,,7] = out[,,7] > 1.358
out0[,,8] = out[,,8] > 1.358
out0[,,9] = out[,,9] > 1.358
out0[,,10] = out[,,10] > 1.358
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
#Plot 1
matplot(delta, 100*power, type = 'b', col = 1:10, lty = c(1,2,3), pch = 1, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c(paste0("m=",1:10)), 
  col = 1:10, lty = c(1,2,3), cex = 0.4)

#Plot 2
matplot(delta, 100*power, type = 'b', col = 1:10, lty = c(1,2,3), pch = 1, main = "Power (Zoom-in version)",
        ylim = c(0, 0.3) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c(paste0("m=",1:10)), 
  col = 1:10, lty = c(1,2,3), cex = 0.4)

###-----------------------------------------------------------
### Conclusion: 
### M = 7 has more power
###-----------------------------------------------------------
