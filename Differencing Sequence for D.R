library(NlcOptim)

m = 1
obj <- function(x){
  m = length(x)-1
  s = 0
  for (g in 1:(m)){
    for (j in g:(m)){
      s = s + x[j+1]^2*x[j-g+1]^2
    }
  }
  s
}

con <- function(x){
  f = NULL
  f = rbind(f, sum(x))  
  f = rbind(f, sum(x^4)-1)
  return(list(ceq = f, c = NULL))
}

x1 = c(0.5, -0.1)
s1 =solnl(x1, objfun = obj, confun= con);s1$par
x2 = c(0.5, -0.1, 0.1)
s2 =solnl(x2, objfun = obj, confun= con);s2$par;obj(s2$par)
x3 = c(0.5, -0.1, 0.1, -0.1)
s3 =solnl(x3, objfun = obj, confun= con);s3$par;obj(s3$par)
x4 = c(0.5, -0.1, 0.1, -0.1, 0.2)
s4 =solnl(x4, objfun = obj, confun= con);s4$par;obj(s4$par)
x5 = c(0.5, -0.1, 0.1, -0.1, 0.1, 0.1)
s5 =solnl(x5, objfun = obj, confun= con);s5$par;obj(s5$par)
x6 = c(0.5, -0.1, 0.1, -0.1, 0.1, 0.1, -0.2)
s6 =solnl(x6, objfun = obj, confun= con);s6$par;obj(s6$par)
x7 = c(0.5, -0.1, 0.1, -0.1, 0.1, 0.1, -0.2, 0.1)
s7 =solnl(x7, objfun = obj, confun= con);s7$par;obj(s7$par)
x8 = c(0.5, -0.1, 0.1, -0.1, 0.1, 0.1, -0.2, 0.1, -0.1)
s8 =solnl(x8, objfun = obj, confun= con);s8$par;obj(s8$par)
x9 = c(1, -0.1, 0.1, -0.1, 0.1, 0.1, -0.2, 0.1, -0.1, 0.2)
s9 =solnl(x9, objfun = obj, confun= con);s9$par;obj(s9$par)
x10 = c(1, -0.1, 0.3, -0.5, 0.1, 0.1, -0.2, 0.1, -0.1, 0.2, 0.1)
s10 =solnl(x10, objfun = obj, confun= con);s10$par;obj(s10$par)


s1$par
s2$par
s3$par
s4$par
s5$par
s6$par
s7$par
s8$par
s9$par
s10$par