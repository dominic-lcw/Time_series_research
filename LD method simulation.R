###------------------------------------
### Simulation
###------------------------------------
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

###------------------------------------------------
# IID (n = 400, k = 1/2, mu = f) m = 1, 3, 5
#out1 (IID, h = 1)

# MA included (AR = 0)
#out 21(ma = -0.5, ar = 0)
#out 22(ma = 0.5, ar = 0)

# AR included (MA = 0)
#out31 (ma = 0, ar = -0.5)
#out32 (ma = 0, ar = 0.5)

# ARMA (MA fixed = -0.5)
#out41 (ma = -0.5, ar = -0.5)
#out42 (ma = -0.5, ar = -0.25)
#out43 (ma = -0.5, ar = 0.25)
#out44 (ma = -0.5, ar = 0.5)

# By-case h
#out51 (ma = -0.5, ar = -0.5)
#out52 (ma = -0.5, ar = 0.25)

###-------------------------------------------------


###---------------------------------------
### Simulation results (1.)
###---------------------------------------
n =400
k = 1/2
mu = f(1:n/n)

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 1, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      z = c(rnorm(n*k),rnorm(n*(1-k))*(1+del))
      x = z+mu
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out1 = out

##---------------------------------------
### Simulation results (2.) MA
###---------------------------------------
n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = 0

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out21 = out

n =400
k = 1/2
mu = f((1:n)/n)
ma1 = 0.5
ar1 = 0

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out22 = out


###---------------------------------------
### Simulation results (3.1-3)
###---------------------------------------
n =400
k = 1/2
mu = f((1:n)/n)
ma1 = 0
ar1 = -0.5

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out31 = out


n =400
k = 1/2
mu = f((1:n)/n)
ma1 = 0
ar1 = 0.5

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out32 = out


###---------------------------------------
### Simulation results (ARMA)
###---------------------------------------
n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = -0.5

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out41 = out

n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = -0.25

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out42 = out

n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = 0.25

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out43 = out

n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = 0.5

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'KS std')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      d1 = diff_seq(x, ds1, 1)
      d2 = diff_seq(x, ds2, 1)
      d3 = diff_seq(x, ds5, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = KS_std(log(d1^2))
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out44 = out


###---------------------------------------
### Best H estimate
###---------------------------------------
n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = -0.5

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'm = 1, h=1')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      h = h_estimate(x)
      d1 = diff_seq(x, ds1, h)
      d2 = diff_seq(x, ds2, h)
      d3 = diff_seq(x, ds5, h)
      d4 = diff_seq(x, ds1, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = dk2(d4^2)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out51 = out

n =400
k = 1/2
mu = f((1:n)/n)
ma1 = -0.5
ar1 = 0.5

if(1){
  n_sim = 2^9
  delta = seq(from = 0, to = 3, length.out = 21)
  out = array(NA, dim =c(n_sim, length(delta), 4),
              dimnames = list(paste0('isim=',1:n_sim),
                              paste0('delta=', delta),
                              c('m = 1', 'm = 3',' m = 5', 'm = 1, h=1')))
  for(i_sim in 1:n_sim){
    set.seed(i_sim)
    for(i_delta in 1:length(delta)){
      del = delta[i_delta]
      a = arma11(n, ma = ma1, ar = ar1)
      z = a*c(rep(1, n*k),rep(3, n*(1-k)))
      x = mu+z
      h = h_estimate(x)
      d1 = diff_seq(x, ds1, h)
      d2 = diff_seq(x, ds2, h)
      d3 = diff_seq(x, ds5, h)
      d4 = diff_seq(x, ds1, 1)
      out[i_sim, i_delta, 1] = dk2(d1^2)
      out[i_sim, i_delta, 2] = dk2(d2^2)
      out[i_sim, i_delta, 3] = dk2(d3^2)
      out[i_sim, i_delta, 4] = dk2(d4^2)
    }
    if(i_sim%%10==0) cat(i_sim, ' >> ')
  }
}
out52 = out


###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###--------------------------------------------------
### Plotting of results
###--------------------------------------------------
###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out0 = out1*0
out0[,,1] = out1[,,1] > 1.358
out0[,,2] = out1[,,2] > 1.358
out0[,,3] = out1[,,3] > 1.358
out0[,,4] = out1[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')

matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "IID", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('m = 1', 'm = 3',' m = 5', 'KS std'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)


out0 = out21*0
out0[,,1] = out21[,,1] > 1.358
out0[,,2] = out21[,,2] > 1.358
out0[,,3] = out21[,,3] > 1.358
out0[,,4] = out21[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')

matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "MA = -0.5", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('m = 1', 'm = 3',' m = 5', 'KS std'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)


out0 = out22*0
out0[,,1] = out22[,,1] > 1.358
out0[,,2] = out22[,,2] > 1.358
out0[,,3] = out22[,,3] > 1.358
out0[,,4] = out22[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')

matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "MA = 0.5", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('m = 1', 'm = 3',' m = 5', 'KS std'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)


###---------------------------------------------
### Part 3
###---------------------------------------------
out0 = out31*0
out0[,,1] = out31[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out31[,,2] > 1.358
out0[,,3] = out31[,,3] > 1.358
out0[,,4] = out31[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')

matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "n=100, k=1/2, mu=f(1:n)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Spline', 'Textbook','log(d^2)','delta KS'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)

out0 = out32*0
out0[,,1] = out32[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out32[,,2] > 1.358
out0[,,3] = out32[,,3] > 1.358
out0[,,4] = out32[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')

matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "n=200, k=1/2, mu=f(1:n)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Spline', 'Textbook','log(d^2)','delta KS'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)

out0 = out33*0
out0[,,1] = out33[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out33[,,2] > 1.358
out0[,,3] = out33[,,3] > 1.358
out0[,,4] = out33[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "n=400, k=1/2, mu=f(1:n)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Spline', 'Textbook','log(d^2)','delta KS'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)

###---------------------------------------------
### Part 4
###---------------------------------------------
out0 = out41*0
out0[,,1] = out41[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out41[,,2] > 1.358
out0[,,3] = out41[,,3] > 1.358
out0[,,4] = out41[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "n=400, k=1/8, mu=f(1:n)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Spline', 'Textbook','log(d^2)','delta KS'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)

out0 = out42*0
out0[,,1] = out42[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out42[,,2] > 1.358
out0[,,3] = out42[,,3] > 1.358
out0[,,4] = out42[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "n=400, k=1/4, mu=f(1:n)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Spline', 'Textbook','log(d^2)','delta KS'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)

out0 = out43*0
out0[,,1] = out43[,,1] > -log(-log(1-0.05)/2)
out0[,,2] = out43[,,2] > 1.358
out0[,,3] = out43[,,3] > 1.358
out0[,,4] = out43[,,4] > 1.358
power = apply(out0, c(2,3), mean)
par(mfrow = c(1, 1))
col = c('blue','red', 'darkgreen','lightgreen')
#Plot 1
matplot(delta, 100*power, type = 'l', col = col, lwd = 2, main = "n=400, k=1/2, mu=f(1:n)", ylim = c(0, 100),
        ylab = expression(K(delta)~"/%"), xlab = expression(delta))
abline(h  = c(0.05, 0, 1)*100, lwd = 0.5, lty = 2)
abline(v = 0, lwd = 0.5, lty = 2)
legend('bottomright', legend = c('Spline', 'Textbook','log(d^2)','delta KS'), 
  col = c('blue','red','darkgreen','lightgreen'), lty = c(1,2,3), cex = 0.6)