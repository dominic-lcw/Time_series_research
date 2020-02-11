f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *10 + (t-0.2)^3 *20
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

