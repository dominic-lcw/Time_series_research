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

ma1 = ma(n, c(0.4,0.5,-0.6,-0.9,0,0,0,0.9), 8)
acf(ma1)
ts.plot(ma1)

n = 10000
#a = ar(1000,0.8, 0.3)
a = ma1
mu = f((1:n)/n)
z = c(rnorm(n/2,1), rnorm(n/2, 1)*5) 
x = z+mu+a
par(mfrow = c(1, 1))
ts.plot(x)

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


d1 = diff_seq(x, m1, 1)
d2 = diff_seq(x, m2, 1)
d3 = diff_seq(x, m3, 1)
d4 = diff_seq(x, m4, 1)
d5 = diff_seq(x, m5, 1)
d6 = diff_seq(x, m6, 1)
d7 = diff_seq(x, m7, 1)
d8 = diff_seq(x, m8, 1)

par(mfrow = c(2, 4))
acf(d1)
acf(d2)
acf(d3)
acf(d4)
acf(d5)
acf(d6)
acf(d7)
acf(d8)

acf(d1^2)
acf(d2^2)
acf(d3^2)
acf(d4^2)
acf(d5^2)
acf(d6^2)
acf(d7^2)
acf(d8^2)

dm1 = log(diff_seq(x, m1, 1)^2)
dm2 = log(diff_seq(x, m2, 1)^2)
dm3 = log(diff_seq(x, m3, 1)^2)
dm4 = log(diff_seq(x, m4, 1)^2)
dm5 = log(diff_seq(x, m5, 1)^2)
dm6 = log(diff_seq(x, m6, 1)^2)
dm7 = log(diff_seq(x, m7, 1)^2)
dm8 = log(diff_seq(x, m8, 1)^2)

acf(dm1)
acf(dm2)
acf(dm3)
acf(dm4)
acf(dm5)
acf(dm6)
acf(dm7)
acf(dm8)

j1 = diff_seq(x, m1, 1)
j2 = diff_seq(x, m1, 2)
j3 = diff_seq(x, m1, 3)
j4 = diff_seq(x, m1, 4)
j5 = diff_seq(x, m1, 5)
j6 = diff_seq(x, m1, 6)
j7 = diff_seq(x, m1, 7)
j8 = diff_seq(x, m1, 8)

acf(j1)
acf(j2)
acf(j3)
acf(j4)
acf(j5)
acf(j6)
acf(j7)
acf(j8)


jm1 = log(diff_seq(x, m1, 1)^2)
jm2 = log(diff_seq(x, m1, 2)^2)
jm3 = log(diff_seq(x, m1, 3)^2)
jm4 = log(diff_seq(x, m1, 4)^2)
jm5 = log(diff_seq(x, m1, 5)^2)
jm6= log(diff_seq(x, m1, 6)^2)
jm7= log(diff_seq(x, m1, 7)^2)
jm8= log(diff_seq(x, m1, 8)^2)

acf(jm1)
acf(jm2)
acf(jm3)
acf(jm4)
acf(jm5)
acf(jm6)
acf(jm7)
acf(jm8)

















