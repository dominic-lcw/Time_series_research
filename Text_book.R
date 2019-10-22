###-----------------------------------------------
### The data. **Only suitable for IID data**
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
z = c(rt(n/2,3), rt(n/2, 3)) 
x = z+mu
par(mfrow = c(1, 1))
ts.plot(x)
d = diff(x) / (sqrt(2))
ts.plot(d)

###-----------------------------------------------
### Calculate test statistics
###-----------------------------------------------
M <- function(x, t){
  n = length(x)
  bar_x = mean(x)
  (sum((x[1:t] - bar_x)^2) - t/n*sum((x - bar_x)^2)) / (sqrt(n))
}

K <- function(x, gam){
  n = length(x)
  k_t = rep(NULL, n)
  for (i in 1:n){
    k_t[i] = M(x, i)/((i/n)*(1-i/n))^gam
  }
  which.min(k_t)
}

V_k <- function(k, x){
  bar = mean(x)
  n = length(x)
  (sum((x-bar)^4) - (1/k)*(sum((x[1:k]-bar)^2)^2) - (1/(n-k))*(sum((x[(k+1):n]-bar)^2, na.rm = TRUE)^2))/n
}

textbook <- function(d){
  k = K(d, 1/3)
  v = V_k(k, d)
  M_t <- rep(NULL, n)
  for (i in 1:(n-1)){
    M_t[i] = M(d, i) / sqrt(v)
  }
  max(abs(M_t))
}

#Testing of code
k = K(d, 1/3)
v = V_k(k, d)
M_t <- rep(NULL, n)
for (i in 1:(n-1)){
  M_t[i] = M(d, i) / sqrt(v)
}
ts.plot(M_t)


###-----------------------------------------------
### Simulation Part
###-----------------------------------------------

if(1){
  n = 500
  n_sim = 80
  delta = seq(from = 0, to = 1, by = 0.05)
  out = array(NA, dim =c(n_sim, length(delta), 1),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('text_book')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      mu = f((1:n/n))
      z = c(rt(n/2,6), rt(n/2, 6)*(1+del))
      x = z+mu
      d = diff(x) / (sqrt(2))
      out[i_sim, i_delta, 1] = textbook(d)
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
power = apply(out0, c(2,3), mean)
power
###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue','red')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend(0.8, 20, legend = c('Textbook','self'), col = c('blue','red'), lty = c(1,2), cex = 0.6)
#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.5) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.2))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)

