rm(list = ls())

###########################################################################
# Estimation Functions
###########################################################################
f.bid_function_fast = function(cost, num_bids, mu, alpha, gamma_1p1oa){
  
  if (exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha) == 0) {
    return(cost + mu/alpha*(num_bids-1)^(-1/alpha)*1/gamma_1p1oa*
             ((num_bids-1)*(gamma_1p1oa/mu*cost)^alpha)^(1/alpha-1))
  } 
  
  cost + 1/alpha*(mu/gamma_1p1oa)*(num_bids-1)^(-1/alpha)*
    pgamma((num_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha)
  # Check gamma(1/alpha) part
}
vf.bid_function_fast = Vectorize(FUN = f.bid_function_fast,vectorize.args = "cost")

vf.w_integrand_z_fast = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, param.u){  
  
  b_z = vf.bid_function_fast(cost=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z
  
  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    dlnorm(u_z, meanlog=(-param.u^2*1/2), sdlog = param.u) # Note: can swap for different distributions
  
  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
  
}

f.funk = function(data_vec, param.u){
  val = integrate(vf.w_integrand_z_fast, w_bid=data_vec[1], 
                  num_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], param.u=param.u, lower=0, upper=Inf,
                  abs.tol = 1e-10)
  #                   rel.tol = .Machine$double.eps^0.7)                    
  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

f.ll_parallel = function(par, y, n, h_x, cl){
  params = par
  v.y = y
  v.n = n
  m.h_x = h_x

  v.mu = params[1]
  v.alpha = params[2]
  param.u = params[3]

  param.h = params[4:(3+dim(m.h_x)[2])]
  v.h = exp(colSums(param.h*t(m.h_x)))
  
  if(param.u <= 0.1) return(-Inf) # Check that these hold at estimated values
  if(sum(v.mu <= 0) > 0) return(-Inf)
  if(sum(v.alpha <= 0.01) > 0) return(-Inf)
  
  # Y Component  
  v.gamma_1p1opa = gamma(1 + 1/v.alpha)
  v.w = v.y/v.h
  dat = cbind(v.w, v.n, v.mu, v.alpha, v.gamma_1p1opa)
  v.f_w = parApply(cl = cl, X=dat, MARGIN=1, FUN=f.funk, param.u=param.u)
  v.f_y = v.f_w/v.h
  return(-sum(log(v.f_y)))
}


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

result
