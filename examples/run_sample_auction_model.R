rm(list = ls())
source("./R/code_sketch.R")


#######################################################
# Load Data
#######################################################

set.seed(301)
# data = # Generate some data
# y, n, x1, x2: positive
# n: discrete and > 1
# y is some function of n, x1, x2

obs = 200
w = rlnorm(obs)
x1 = rlnorm(obs) + .5*w
x2 = .1*rlnorm(obs) + .3*w
e = 2*rlnorm(obs)
n = sample(2:10, obs, replace=TRUE)
y = 10 - .5*n + x1 + x2 + e
data = data.frame(cbind(y, n, x1, x2))
plot(n, y)

v.y = data$y
v.n = data$n
m.h_x = as.matrix(cbind(log(data$x1),log(data$x2)))

# inital parameter guess
x0 =  c(8, 2, .5, .4, .6)

library(parallel)
cl = makeCluster(4)
clusterExport(cl,varlist=c("vf.bid_function_fast",
                           "vf.w_integrand_z_fast",
                           "f.funk"))

f.ll_parallel(x0, y = v.y, n = v.n, h_x = m.h_x, cl = cl)

optim_control = list(maxit = 2000, parscale = c(1, 0.1, 1, 0.1,
                                                rep(1, length(x0) - 4)))

result = optim(par = x0, fn = f.ll_parallel, control=optim_control,
               y=v.y, n=v.n, h_x=m.h_x, cl=cl)

stopCluster(cl = cl)

result
