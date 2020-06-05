#rm(list = ls())
#source("./R/code_sketch.R")

library(auctionr)

#######################################################
# Load Data
#######################################################

set.seed(301)
# data = # Generate some data
# y, n, x1, x2: positive
# n: discrete and > 1
# y is some function of n, x1, x2

obs = 50
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

#### Testing the package on data generated using a wrong model

res <- auction_model(data,
                     #init_param =  x0,
                     init_param =  c(12, 4.7, 0.2, 0.07, 0.05),
                     num_cores = 25,
                    # method = "Nelder-Mead",
                     control = list(maxit = 2000),
                     std_err = TRUE)

### Temporary code to insert

library(maxLik)

res_ML <- maxLik(start=init_param, logLik=f__ll_parallelneg,
      y=v__y, n=v__n, h_x=m__h_x, cl=cl) #method=

print(res_ML)
print(stdEr(res_ML))


### Default Example

set.seed(100)
dat <- auction_generate_data(obs = 100, mu = 10, alpha = 2, sigma = 0.2,
                             beta = c(-1,1), new_x_mean= c(-1,1), new_x_sd = c(0.5,0.8))

res <- auction_model(dat,
                     init_param =  c(8, 2, .5, .4, .6),
                     num_cores = 1,
                     method = "BFGS",
                     control = list(trace=1, parscale = c(1,0.1,0.1,1,1)),
                     std_err = TRUE)

## Example previously not working

set.seed(100)

dat <- auction_generate_data(obs = 100, mu = 10, alpha = 2, sigma = 0.2,
                             beta = c(-1,1), new_x_mean= c(-1,1), new_x_sd = c(0.5,0.8))

auction_model(dat,
              init_param =  c(1, 1, 1, 0, 0),
              num_cores = 1,
              method = "BFGS",
              control = list(trace=1, parscale = c(1,0.1,0.1,1,1)),
              std_err = TRUE)
