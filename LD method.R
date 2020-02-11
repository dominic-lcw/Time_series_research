###-----------------------------------------------
### The data. **AR included**
###-----------------------------------------------
f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *10 + (t-0.2)^3 *20
}
arma11 = function(n, ma, ar){
  wn = rnorm(n)
  arma11  = rep(NA, n)
  arma11[1] = wn[1]
  for (i in 2:n){
    arma11[i] = arma11[i-1]*ar + wn[i-1]*ma + wn[i]
  }
  arma11
}
n = 400
mu = f((1:n)/n)
a = arma11(n, ma = -0.8, ar = -0.4)
z = a*c(rep(1, n/2),rep(3, n/2))
x = mu+z
ts.plot(x)

###----------------------------------
### LD2 Method
###----------------------------------
ds1 = s1$par
ds2 = s2$par
ds3 = s3$par
ds4 = s4$par
ds5 = s5$par
ds6 = s6$par
ds7 = s7$par
ds8 = s8$par
ds9 = s9$par
ds10 = s10$par


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

h1 = c(0.7071, -0.7071)
h2 = c(0.8090, -0.5, -0.3090)
h3 = c(0.1942, 0.2809, 0.3832, -0.8582)

h_estimate = function(x){
  d = diff_seq(x, h1, 1)
  autocorr = acf(d, plot = FALSE)$acf
  tr = which((autocorr < sqrt(1.96/100)) ==TRUE)
  if (length(tr) !=0){
    h = min(tr)
  }else{
    h = which(autocorr == min(autocorr))
  }
  h
}


kern = function(x){ #Bartlett kernel
  ifelse(abs(x)<=1,1-abs(x), 0)
}


#This method fails
dk2.lrv = function(x, l){
	n = length(x)
	l = floor(l)
	sum = 0
	covar = c(acf(x, l+1, type = 'cov', plot=FALSE)$acf)
	k = (-l):l
	out = sum(kern(k/l) * covar[abs(k) +1]) /mean(x)^2
	return(out)
}


dk2 = function(x){ #x is the d^2, dseq is the differencing sequence
  n = length(x)
  s0 = rep(NA, n-1)
  c = dk2.lrv(x, 2*floor(n^(1/3)))

  for (i in 1:(n-1)){
	  m1 = log(mean(x[1:i]))
	  m2 = log(mean(x[(i+1):n]))
	  s0[i] = (m1-m2)*(i/n)*(1-i/n)
  } 
  max(abs(s0))*sqrt(n)/sqrt(c)
}
ts.plot(d2/mean(d2))

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

ts.plot(a*c(rep(1, n/2),rep(1+del, n/2)))
z = c(rep(1, n/2),rep(1+del, n/2))

if(1){
  n = 400
  n_sim = 2^6
  delta = seq(from = 0, to = 2, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta),3),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('delta h1','delta h3','delta h8')))
  mu = f((1:n)/n)
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      z = arma11(n, -0.8, -0.4)*c(rep(1, n/2),rep(1+del, n/2))
      # z = c(rnorm(n/2,0,1),rnorm(n/2,0,1)*(1+del))
      x = z+mu#+a
      h = h_estimate(x)
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds1,h)
      d3 = diff_seq(x, ds5,h)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
    }
    if(i_sim%%50==0) cat(i_sim, ' >> ')
  }
}

par(mfrow = c(1, 2))
out0 = out*0
out0[,,1] = out[,,1] > 1.358
out0[,,2] = out[,,2] > 1.358
out0[,,3] = out[,,3] > 1.358
power = apply(out0, c(2,3), mean)
power
#Plot 1
matplot(delta, 100*power, type = 'b', col = 1:10, lty = c(1,2,3), 
  pch = "12359K", main = "Power (n = 400, ARMA(1,1))", ylim = c(0, 100),
    ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('h=1','h=3', 'h by case'), 
  col = 1:10, lty = c(1,2,3), cex = 0.8)

matplot(delta, 100*power,  type = 'b', col = 1:10, lty = c(1,2,3), 
  pch = "12359K", main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)


###---------------------------------------
### Multiple sequence
###---------------------------------------

if(1){
  n = 400
  n_sim = 2^10
  delta = seq(from = 0, to = 1, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 6),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('delta 1','delta 2','delta 3','delta 5','delta 9','KS')))
  mu = f((1:n)/n)
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      # z = arma11(n, -0.8, -0.4)*c(rep(1, n/2),rep(1+del, n/2))
      z = c(rnorm(n/2,0,1),rnorm(n/2,0,1)*(1+del))
      x = z+mu#+a
      out[i_sim, i_delta, 1] = dk2(diff_seq(x, s1$par,1)^2)
      out[i_sim, i_delta, 2] = dk2(diff_seq(x, s2$par,1)^2)
      out[i_sim, i_delta, 3] = dk2(diff_seq(x, s3$par,1)^2)
      out[i_sim, i_delta, 4] = dk2(diff_seq(x, s5$par,1)^2)
      out[i_sim, i_delta, 5] = dk2(diff_seq(x, s9$par,1)^2)
      out[i_sim, i_delta, 6] = KS_std(log((diff(x)/sqrt(2))^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}

out0 = out*0
out0[,,1] = out[,,1] > 1.358
out0[,,2] = out[,,2] >1.358
out0[,,3] = out[,,3] >1.358
out0[,,4] = out[,,4] >1.358
out0[,,5] = out[,,5] >1.358
out0[,,6] = out[,,6] >1.358
power = apply(out0, c(2,3), mean)
power

par(mfrow = c(1, 2))
#Plot 1
matplot(delta, 100*power, type = 'b', col = 1:10, lty = c(1,2,3), 
	pch = "12359K", main = "Power", ylim = c(0, 100),
    ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('m=1','m=2','m=3','m=5','m=9','KS'), 
  col = 1:10, lty = c(1,2,3), cex = 0.5)

matplot(delta, 100*power,  type = 'b', col = 1:10, lty = c(1,2,3), 
	pch = "12359K", main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)


