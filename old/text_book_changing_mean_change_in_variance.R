###-----------------------------------------------
### Change in variance in changing mean
###-----------------------------------------------
f = function(t){
  5*cos(t*9) - (t - 0.7)^2 *20 + (t-0.2)^3 *10
}
n = 1000
mu = f((1:n)/n)
z = c(rnorm(n/2)*1, rnorm(n/2)*4)
x = z+mu

#Plot the sample graph
par(mfrow = c(1,2))
ts.plot(mu)
ts.plot(x)

###-----------------------------------------------
### Change in variance in changing mean
###-----------------------------------------------

bar_x = mean(x)
sq = (x - bar_x)^2
bar_sq = mean(sq)
Tn = cumsum(sq - bar_sq) / sqrt(n)
ts.plot(abs(Tn) / sd(sq)) #Using sd(.) as the proxy for sd


###-----------------------------------------------------------
### Simulations
###-----------------------------------------------------------

if(1){
  n = 1000
  n_sim = 50
  delta = seq(from = 0, to = 2, by = 0.05)
  out = array(NA, dim =c(n_sim, length(delta), 2),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('std', 'sn')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      d = delta[i_delta]
      mu = f((1:n/n))
      z = c(rnorm(n/2)*1, rnorm(n/2) * (1+d))
      x = z+mu
      d1 = diff(x)
      sq = (d1 - mean(d1))^2
      out[i_sim, i_delta, 1] = KS_std(sq)
      out[i_sim, i_delta, 2] = KS_sn(sq)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}


###-----------------------------------------------------------
### Power
###-----------------------------------------------------------
out0 = out*0
out0[,,1] = out[,,1] > 1.358
out0[,,2] = out[,,2] > 40.1
power = apply(out0, c(2,3), mean)
power

###-----------------------------------------------------------
### Plot results
###-----------------------------------------------------------
par(mfrow = c(1, 2))
col = c('blue', 'red')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
#Plot 2
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "Power (Zoom-in version)",
        ylim = c(0, 0.2) *100, ylab = expression(K(delta)~"/%"), xlab = expression(delta),
        xlim = c(0, 0.5))
abline(h = c(0:5), lwd = 0.5, lty = 3)
abline(v = 0, lwd = 0.5, lty = 2)
















