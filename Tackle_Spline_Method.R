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
### Testing Cases
###-----------------------------------------------
n = 1000
del = 1
z = c(rnorm(n/2,1),rnorm(n/2,1)*(1+del)) #Normal Distributio
z = c(rt(n/2,6), rt(n/2, 6)*(1+del)) #T-distribution
z = c(rexp(n/2, 1), rexp(n/2, 1)*(1+del))

#Bimodal Distribution
p1 = c(rbinom(n*1/2, 1, 0.5))
p2 = c(rbinom(n*1/2, 1, 0.5))
z = c(rnorm(n*1/2,-10)*p1+rnorm(n*1/2,10)*(1-p1),
	(rnorm(n*1/2,-10, (1+del))*p2+rnorm(n*1/2,10,(1+del))*(1-p2))) #Bimodal T-distribution

z= c(rnorm(n/2,1),rnorm(n/4,1,(2+del)), rnorm(n/4,1,(4+del))) #Normal Distribution
par(mfrow = c(1,1))
ts.plot(z)

par(mfrow = c(1,2))
hist(z[1:n/2], breaks = 100)
hist(z[n/2+1:n], breaks = 100)

###-----------------------------------------------
### Testing Normal Case
###-----------------------------------------------

###Input Parameter###
n = 400
n_sim = 200

if(1){
 	delta = seq(from = 0, to = 1, length.out = 11)
 	out = array(NA, dim =c(n_sim, length(delta), 5),
                dimnames = list(paste0('isim=',1:n_sim),
                paste0('delta=', delta),c('Normal', 'Student-T', 'Exp', 'Bimodal', '2cp_norm')))
    for(i_sim in 1:n_sim){
        set.seed(i_sim)
        for(i_delta in 1:length(delta)){
            del = delta[i_delta]
      		mu = 0
			#Random Number
			z_norm = c(rnorm(n/2,1),rnorm(n/2,1,(1+del))) #Normal Distribution
			z_t = c(rt(n/2,6), rt(n/2, 6)*(1+del))# T-distribution
			z_exp = c(rexp(n/2, 1)*1, rexp(n/2, 1)*(1+del)) #Exponential Distribution
			z_bi = c(rnorm(n*1/2,-2,1)*p1+rnorm(n*1/2,10,1)*(1-p1),
				(rnorm(n*1/2,-10,(1+del))*p2+rnorm(n*1/2,10,(1+del))*(1-p2))) #Bimodal T-distribution
			z_norm = c(rnorm(n/2,1),rnorm(n/4,1)*(1+del), rnorm(n/4,1)*(4+del)) #Normal Distribution
			z_2cp_norm= c(rnorm(n/2,1),rnorm(n/4,1,(2+del)), rnorm(n/4,1,(4+del))) #Multiple Change Point
			
			# #Time Series
			x_norm = z_norm+mu
			x_t = z_t + mu
			x_exp = z_exp + mu
			x_bi = z_bi + mu
			x_2cp_norm = z_2cp_norm + mu
			
			# #Test Statistics
			out[i_sim, i_delta, 1] = spline_test(x_norm)
			out[i_sim, i_delta, 2] = spline_test(x_t)
			out[i_sim, i_delta, 3] = spline_test(x_exp)
			out[i_sim, i_delta, 4] = spline_test(x_bi)
			out[i_sim, i_delta, 5] = spline_test(x_2cp_norm)
		}
    	if(i_sim%%10==0) cat(i_sim, ' >> ')
  	}
}

###-----------------------------------------------
### Result
###-----------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out[,,2] > -log(-log(1-0.05)/2)
out0[,,3] = out[,,3] > -log(-log(1-0.05)/2)
out0[,,4] = out[,,4] > -log(-log(1-0.05)/2)
out0[,,5] = out[,,5] > -log(-log(1-0.05)/2)
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
#Plot 1
matplot(delta, 100*power, type = 'b', col = 1:5, lty = c(1,2,3), pch = 1, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Normal', 'Student-T', 'Exp', 'Bimodal', 'Multi_cp'), 
  col = 1:5, lty = c(1,2,3), cex = 0.7)

#Plot 2
matplot(delta, 100*power, type = 'b', col = 1:5, lty = c(1,2,3), pch = 1, main = "Power (Zoom-in version)",
        ylim = c(0, 0.3) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Normal', 'Student-T', 'Exp', 'Bimodal', 'Multi_cp'), 
  col = 1:5, lty = c(1,2,3), cex = 0.4)

