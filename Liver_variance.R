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
d = diff(x) / (sqrt(2))
ts.plot(d)

###-----------------------------------------------
### The data. **AR included**
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
  if (T> -log(-log(1-alpha))/2){
      return(TRUE)
  }else{
    return(FALSE)
  }
}

spline_test = function(x){
  epsilon = 0.00001
  n = length(x)
  ax = (1:n)/n
  old_d = 0
  count = 0
  while (count<50){
    g = smooth.spline(ax, x, w = weight, cv = TRUE)$y
    result= MaxL(g, fit$y)
    new_d = result[1]
    tau = result[2]
    if ((new_d - old_d)<epsilon){
      return(Liver_test(new_d, 0.05))
    }else{
      sigma1 = MLsigma1(x, fit$y, tau)
      sigma2 = MLsigma2(x, fit$y, tau)
      weight = c(rep((1/sigma1)^2, n/2), rep((1/sigma2)^2, n/2))
    }
    count = count+1
    old_d = new_d  
  }
  return(FALSE)
}

