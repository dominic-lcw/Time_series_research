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
ts.plot(d2)
mean(d2[1:200])
mean(d2[201:399])

ts.plot(artificial_x(d2))
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

KS_std_new = function(d2){ #
  k = k_split(d2)
  sd = two_sigma(d2,k)
  n = length(d2)
  d2_bar = mean(d2)
  Tn =cumsum(d2 - d2_bar) /sqrt(c(rep(sd[1],k),rep(sd[2],(n-k))))
  TS = max(abs(Tn) / sqrt(n))
)
k
sd

ts.plot(d2)
abline(h = d2_bar, col = 'red', lwd = 1.5)


ts.plot(Tn)
var_head(d2, 350)

ts.plot(x)
a= KS_std_new(d2)
b = KS_std(d2)
a
b
ts.plot(d2)
ts.plot(Tn)
x
rep(0,3.5)

  ?ceiling