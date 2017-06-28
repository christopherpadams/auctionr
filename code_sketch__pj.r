library(parallel)

auction <- function() {
  
  # Check inputs
  
  
  # Prepare data
  dat = genSampleData1()
  
  v__y = dat$v__y
  v__n = dat$v__n
  m__h_x = dat$m__h_x
   
  # Get initial guess for convergence
  x0 = f_conv_initGuess(v__y, v__n, m__h_x)

  # Set up parallelization
  cl = makeCluster(10)
  clusterExport(cl,varlist=c("vf__bid_function_fast",
                             "vf__w_integrand_z_fast",
                             "f__funk"))
  # f__ll_parallel(x0, y = v__y, n = v__n, h_x = m__h_x, cl = cl)

  #   Get convergence conditions
  optim_control = f_conv_conditions(x0)
  
  # Run
  result = optim(par=x0, fn=f__ll_parallel, control=optim_control,
                 y=v__y, n=v__n, h_x=m__h_x, cl=cl)
  
  # Inspect result
  
  # Return result
  return(result)
  return("test")
}



genSampleData1 <- function() {
  # For testing purposes, we will generate sample data
  
  set.seed(301)
  # data = # Generate some data
  # y, n, x1, x2: positive
  # n: discrete and > 1
  # y is some function of n, x1, x2
  
  obs = 200
  w = rlnorm(obs)
  x1 = rlnorm(obs) + 0.5*w
  x2 = 0.1*rlnorm(obs) + 0.3*w
  e = 2*rlnorm(obs)
  n = sample(2:10, obs, replace=TRUE)  
  y = 10 - 0.5*n + x1 + x2 + e
  data = data.frame(cbind(y, n, x1, x2))
  #plot(n, y)
  
  v__y = data$y
  v__n = data$n
  m__h_x = as.matrix(cbind(log(data$x1),log(data$x2)))
  
  return(
    list(
      v__y = v__y, v__n = v__n, m__h_x = m__h_x
      )
  )
}



f_conv_initGuess <- function(v__y, v__n, m__h_x) {
  # Find starting point / intial guess
  #
  # Input:
  #   v__y
  #
  #   v__n
  #
  #   m__h_x
  #
  # Output:
  #   Returned list is of the form
  #     (mu, a, u, h[1], h[2], ... )
  
  # Currently, fix initial guess
  mu_init = 8
  a_init = 2
  u_init = 0.5
  h_init = rep(0.5, dim(m__h_x)[2])
  h_init[1] = 0.4
  h_init[2] = 0.5
  
  return(c(mu_init, a_init, u_init, h_init))
}



f_conv_conditions <- function(x_initialGuess) {
  # Define convergence conditions
  #
  # Input:
  #   x_initialGuess
  #     must be a numeric array of at least length 5
  return(
    list(
      maxit = 2000,
      parscale = c(1, 0.1, 1, 0.1, rep(1, length(x_initialGuess) - 4) ))
    )
}



###########################################################################
# Estimation Functions
###########################################################################
f__bid_function_fast = function(cost, num_bids, mu, alpha, gamma_1p1oa){

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
vf__bid_function_fast = Vectorize(FUN = f.bid_function_fast,vectorize.args = "cost")

vf__w_integrand_z_fast = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, param.u){

  b_z = vf__bid_function_fast(cost=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z

  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    dlnorm(u_z, meanlog=(-param.u^2*1/2), sdlog = param.u) # Note: can swap for different distributions

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)

}



f__funk = function(data_vec, param_u){
  
  val = integrate(vf__w_integrand_z_fast, w_bid=data_vec[1],
                  num_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], param_u=param_u, lower=0, upper=Inf,
                  abs.tol = 1e-10)
  #                   rel.tol = .Machine$double.eps^0.7)
  if(val$message != "OK")
    stop("Integration failed.")
  return(val$value)
}



f__ll_parallel = function(x0, y, n, h_x, cl){
  params = x0
  v__y = y
  v__n = n
  m__h_x = h_x

  v__mu = params[1]
  v__alpha = params[2]
  u = params[3]

  h = params[4:( 3 + dim(m__h_x)[2] )]
  v__h = exp( colSums( h * t(m__h_x) ) )

  if (u <= 0.1)
    return(-Inf) # Check that these hold at estimated values
  else if ( sum (v__mu <= 0 ) > 0 )
    return(-Inf)
  else if ( sum( v.alpha <= 0.01 ) > 0)
    return(-Inf)
  else
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/v__alpha)
    v__w = v__y / v__h
    dat = cbind(v__w, v__n, v__mu, v__alpha, v__gamma_1p1opa)
    
    v__f_w = parApply(cl = cl, X = dat, MARGIN = 1, FUN = f__funk, param_u = u)
    v__f_y = v__f_w / v__h
    return(-sum(log(v.f_y)))
}






###########################################################################
# Run code
###########################################################################
res = auction()
print(res)