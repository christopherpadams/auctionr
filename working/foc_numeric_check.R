hEnv_tmp = new.env()
x0 =  c(7.93369419885821, 0.66080453704268, 0.748151129180016, 0.399372886551224, -0.00062711344877579)
dat_X = dat_X=t( data[ ! names(data) %in% c("winning_bid", "n_bids") ] )
dat__winning_bid = data$winning_bid
dat__n_bids = data$n_bids
cl = parallel::makeCluster(2)
listFuncCall = list(funcName = "dgamma",
                    funcID = auctionmodel:::auction__get_id_distrib(sFuncName = "dgamma"))

parallel::clusterExport(cl,
                        varlist=c("vf__bid_function_fast",
                                  "vf__w_integrand_z_fast",
                                  "f__funk",
                                  "auction__deriv_bid_alpha_integrand",
                                  "auction__deriv_bid_alpha_integrate",
                                  "auction__deriv_bid_alpha",
                                  "auction__vec_deriv_bid_alpha",
                                  "auction__deriv_pv_alpha",
                                  "auction__deriv_alpha_integrand",
                                  "auction__deriv_alpha_integrate",
                                  "auction__deriv_bid_mu_integrand",
                                  "auction__deriv_bid_mu_integrate",
                                  "auction__deriv_bid_mu",
                                  "auction__vec_deriv_bid_mu",
                                  "auction__deriv_pv_mu",
                                  "auction__deriv_mu_integrand",
                                  "auction__deriv_mu_integrate",
                                  "auction__deriv_sig_integrate",
                                  "auction__deriv_sig_integrand",
                                  "auction__deriv_beta_integrand",
                                  "auction__deriv_beta_integrate",
                                  "auction__deriv_ll_beta",
                                  "auction__deriv_gradient",
                                  "auction__deriv_unobs",
                                  "auction__f_unobs_gamma",
                                  "auction__f_unobs_lognorm",
                                  "auction__f_unobs_weibull"),
                        envir = environment(auction_model) )


grad(func = f__ll_parallel,
     x = x0,
     method = "Richardson",
     dat__winning_bid = dat__winning_bid,
     dat__n_bids = dat__n_bids,
     dat_X = dat_X,
     listFuncCall = listFuncCall,
     hTracker = hTracker,
     cl = cl)


auction__deriv_gradient(x0 = x0,
                        dat__winning_bid = dat__winning_bid,
                        dat__n_bids = dat__n_bids,
                        dat_X = dat_X,
                        cl = cl,
                        listFuncCall = listFuncCall,
                        hTracker = hTracker)

parallel::stopCluster(cl)

######## INDIVIDUAL COMPONENTS TESTS
##### UH
####Gamma
###sigma

distrib_std_dev = 0.6
w_bid = 8
n_bids = 5
mu = 5
alpha = 3
gamma_1p1oa = gamma(1+1/alpha)
cost = auction__generate_cost(obs = 300,
                              v.n = rep(5, obs),
                              mu = mu,
                              alpha = alpha)[2]
b_z = vf__bid_function_fast(cost=cost, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)

auction__f_unobs_gamma_sig <- function(distrib_std_dev, w_bid, b_z){

  deriv_unobs_sig = -(1/distrib_std_dev^2)^(2+(1/distrib_std_dev^2))*1/(gamma(1/distrib_std_dev^2))*(w_bid/b_z)^(1/distrib_std_dev^2-1)*
    exp(-(1/distrib_std_dev^2)*w_bid/b_z)*w_bid/b_z*(
      log(1/distrib_std_dev^2) + 1 - digamma(1/distrib_std_dev^2) + log(w_bid/b_z) - w_bid/b_z
    )

  deriv_unobs_u = -1/gamma(1/distrib_std_dev^2)*(1/distrib_std_dev^2)^(1+1/distrib_std_dev^2)*
    (w_bid/b_z)^(1/distrib_std_dev^2 -2)*(distrib_std_dev^2-1+w_bid/b_z)*
    exp(-1/distrib_std_dev^2 * w_bid/b_z)

  f_unobs = (1/(distrib_std_dev^2))^(1/(distrib_std_dev^2))/gamma(1/(distrib_std_dev^2))*
    (w_bid/b_z)^(1/(distrib_std_dev^2)-1)*exp(-1/(distrib_std_dev^2)*w_bid/b_z)
  return(deriv_unobs_sig)
}

auction__f_unobs_gamma_u <- function(distrib_std_dev, w_bid, b_z){

  deriv_unobs_sig = -(1/distrib_std_dev^2)^(2+(1/distrib_std_dev^2))*1/(gamma(1/distrib_std_dev^2))*(w_bid/b_z)^(1/distrib_std_dev^2-1)*
    exp(-(1/distrib_std_dev^2)*w_bid/b_z)*w_bid/b_z*(
      log(1/distrib_std_dev^2) + 1 - digamma(1/distrib_std_dev^2) + log(w_bid/b_z) - w_bid/b_z
    )

  deriv_unobs_u = -1/gamma(1/distrib_std_dev^2)*(1/distrib_std_dev^2)^(1+1/distrib_std_dev^2)*
    (w_bid/b_z)^(1/distrib_std_dev^2 -2)*(distrib_std_dev^2-1+w_bid/b_z)*
    exp(-1/distrib_std_dev^2 * w_bid/b_z)

  f_unobs = (1/(distrib_std_dev^2))^(1/(distrib_std_dev^2))/gamma(1/(distrib_std_dev^2))*
    (w_bid/b_z)^(1/(distrib_std_dev^2)-1)*exp(-1/(distrib_std_dev^2)*w_bid/b_z)

  return(deriv_unobs_u)
}

auction__f_unobs_gamma_f <- function(distrib_std_dev, x){

  f_unobs = (1/(distrib_std_dev^2))^(1/(distrib_std_dev^2))/gamma(1/(distrib_std_dev^2))*
    (x)^(1/(distrib_std_dev^2)-1)*exp(-1/(distrib_std_dev^2)*x)

  return(f_unobs)
}

grad(func = auction__f_unobs_gamma_f,
     x = distrib_std_dev,
     method = "simple",
     w_bid = w_bid,
     b_z = b_z)


auction__f_unobs_gamma_sig(distrib_std = distrib_std_dev,
                           w_bid = w_bid,
                           b_z = b_z)


###first term (w_bid/b_z) - checks out
grad(func = auction__f_unobs_gamma_f,
     x = w_bid/b_z,
     method = "simple",
     distrib_std_dev = distrib_std_dev)


auction__f_unobs_gamma_u(distrib_std = distrib_std_dev,
                           w_bid = w_bid,
                           b_z = b_z)


####Log normal
###sigma^2

auction__f_unobs_lognorm_sig <- function(distrib_std_dev, w_bid, b_z){

  deriv_unobs_sig = -distrib_std_dev^2*1/(w_bid/b_z*sqrt(2*pi))*1/(2*(log(1+distrib_std_dev^2))^(3/2))*1/(1+distrib_std_dev^2) *
    exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))
  -1/(w_bid/b_z*sqrt(2*pi*log(1+distrib_std_dev^2)))*exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))*
    (
      2*(log(1+distrib_std_dev^2)-1)*(log(w_bid/b_z)-1/2*log(1+distrib_std_dev^2))/(4*(log(1+distrib_std_dev^2))^2*(1+distrib_std_dev^2))
    )

  return(deriv_unobs_sig)
}

auction__f_unobs_lognorm_f <- function(distrib_std_dev, w_bid, b_z){

  f_unobs = 1/(w_bid/b_z*sqrt(2*pi*log(1+distrib_std_dev^2)))*
    exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))
  return(f_unobs)
}

grad(func = auction__f_unobs_lognorm_f,
     x = distrib_std_dev,
     method = "simple",
     w_bid = w_bid,
     b_z = b_z)


auction__f_unobs_lognorm_sig(distrib_std = distrib_std_dev,
                           w_bid = w_bid,
                           b_z = b_z)


###first term - checks out
auction__f_unobs_lognorm_f <- function(distrib_std_dev, x){

  f_unobs = 1/(x*sqrt(2*pi*log(1+distrib_std_dev^2)))*
    exp(-(log(x)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))
  return(f_unobs)
}


auction__f_unobs_lognorm_u <- function(distrib_std_dev, w_bid, b_z){

  deriv_unobs_u = -exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))*
    (2*log(w_bid/b_z)+3*log(1+distrib_std_dev^2))/(2*(w_bid/b_z)^2*(log(1+distrib_std_dev^2))^(3/2)*sqrt(2*pi))

  return(deriv_unobs_u)
}

grad(func = auction__f_unobs_lognorm_f,
     x = w_bid/b_z,
     method = "simple",
     distrib_std_dev = distrib_std_dev)

auction__f_unobs_lognorm_u(distrib_std = distrib_std_dev,
                             w_bid = w_bid,
                             b_z = b_z)


####Weibull
###Sigma^2
auction__f_unobs_weibull_sig <- function(distrib_std_dev, w_bid, b_z){

  listParam = auction__get_unobs_params(distrib_std_dev, 3)
  shape = listParam$shape

  dshape = shape^2*gamma(1+1/shape)/(
    digamma(1+1/shape)*gamma(1+2/shape)-2*digamma(1+2/shape)*gamma(1+2/shape)
  )

  df = gamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  - 1/shape*gamma(1+1/shape)*digamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  + shape*gamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(shape - 1)*(
    log(w_bid/b_z*gamma(1+1/shape)) - (shape - 1)*digamma(1+1/shape)/(shape^2)
  )*exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  - shape*gamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(2*shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)*(
      log(w_bid/b_z*gamma(1+1/shape)) - digamma(1+1/shape)/shape
    )

  deriv_unobs_sig =df*dshape

  return(deriv_unobs_sig)
}

auction__f_unobs_weibull_f <- function(distrib_std_dev, w_bid, b_z){

  listParam = auction__get_unobs_params(distrib_std_dev, 3)
  shape = listParam$shape


  f_unobs = shape*gamma(1+1/shape)*(w_bid/b_z*gamma(1+1/shape))^(shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  return(f_unobs)
}


grad(func = auction__f_unobs_weibull_f,
     x = distrib_std_dev,
     method = "simple",
     w_bid = w_bid,
     b_z = b_z)

auction__f_unobs_weibull_sig(distrib_std_dev = distrib_std_dev,
                             w_bid = w_bid,
                             b_z = b_z)

###First term - checks out
auction__f_unobs_weibull_u <- function(distrib_std_dev, w_bid, b_z){
  listParam = auction__get_unobs_params(distrib_std_dev, 3)
  shape = listParam$shape

  deriv_unobs_u = shape*(gamma(1+1/shape))^shape*exp(-(w_bid/b_z*gamma(1+1/shape))^shape)*(
    (shape - 1)*(w_bid/b_z)^(shape -2) - (w_bid/b_z)^(2*shape-2)*shape*(gamma(1+1/shape))^shape
  )

  return(deriv_unobs_u)
}


auction__f_unobs_weibull_f <- function(distrib_std_dev, x){
  listParam = auction__get_unobs_params(distrib_std_dev, 3)
  shape = listParam$shape

  f_unobs = shape*gamma(1+1/shape)*(x*gamma(1+1/shape))^(shape - 1)*
    exp(-(x*gamma(1+1/shape))^shape)
  return(f_unobs)
}


grad(func = auction__f_unobs_weibull_f,
     x = w_bid/b_z,
     method = "simple",
     distrib_std_dev = distrib_std_dev)

auction__f_unobs_weibull_u(distrib_std_dev = distrib_std_dev,
                             w_bid = w_bid,
                             b_z = b_z)


##### ALPHA
distrib_std_dev = 0.6
w_bid = 8
n_bids = 5
mu = 5
alpha = 3
gamma_1p1oa = gamma(1+1/alpha)

###Bid function - checks out
auction__deriv_bid_alpha_integrand <- function(n_bids, gamma_1p1oa, alpha, mu, xi){
  val = -exp(-(n_bids - 1)*(gamma_1p1oa/mu*xi)^alpha)*(n_bids - 1)*(gamma_1p1oa/mu*xi)^alpha*(
    log(gamma_1p1oa/mu*xi) - digamma(1+1/alpha)/alpha
  )
  val[(gamma_1p1oa/mu)^alpha == Inf] = 0
  val[exp(-(n_bids-1)*(gamma_1p1oa/mu*xi)^alpha) == 0] = 0
  return(val)
}

auction__deriv_bid_alpha_integrate <- function(n_bids, gamma_1p1oa, alpha, mu, z){
  int = rep(NA, length(z))
  for (i in 1: length(z)){
    int[i]  = stats::integrate(auction__deriv_bid_alpha_integrand,
                               n_bids = n_bids,
                               gamma_1p1oa = gamma_1p1oa,
                               alpha = alpha,
                               mu = mu,
                               lower = z[i],
                               upper = Inf,
                               abs.tol = 1e-10)$value
  }
  return(int)
}

auction__deriv_bid_alpha <- function(alpha, mu, gamma_1p1oa, z, n_bids){

  alpha_bid_1 = 1/alpha*mu/gamma_1p1oa*(n_bids - 1)^(-1/alpha)*pgamma((n_bids- 1)*(gamma_1p1oa/mu*z)^alpha, 1/alpha, lower = FALSE)*gamma(1/alpha)
  alpha_bid_2 = exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)


  deriv_alpha_bid_1 = auction__deriv_bid_alpha_integrate(n_bids = n_bids,
                                                         gamma_1p1oa = gamma_1p1oa,
                                                         alpha = alpha,
                                                         mu = mu,
                                                         z = z)
  deriv_alpha_bid_2 = exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)*
    (n_bids - 1)*(gamma_1p1oa/mu*z)^alpha*(
      log(gamma_1p1oa/mu*z) - digamma(1+1/alpha)/alpha
    )

  deriv_bid_alpha = deriv_alpha_bid_1*alpha_bid_2 +
    alpha_bid_1*deriv_alpha_bid_2


  return(deriv_bid_alpha)
}

f__bid_function_fast = function(cost, n_bids, mu, alpha){

  if (exp(-(n_bids-1)*(1/(mu/gamma(1+alpha))*cost)^alpha) == 0) {
    return(cost + mu/alpha*(n_bids-1)^(-1/alpha)*1/gamma(1+1/alpha)*
             ((n_bids-1)*(gamma(1+1/alpha)/mu*cost)^alpha)^(1/alpha-1))
  }

  cost + 1/alpha*(mu/gamma(1+1/alpha))*(n_bids-1)^(-1/alpha)*
    stats::pgamma((n_bids-1)*(1/(mu/gamma(1+1/alpha))*cost)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(n_bids-1)*(1/(mu/gamma(1+1/alpha))*cost)^alpha)
}

grad(f__bid_function_fast,
     x = alpha,
     cost = 2,
     n_bids = n_bids,
     mu = mu)
auction__deriv_bid_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z = 2, n_bids = n_bids)


###Private values - checks out

auction__deriv_pv_alpha <- function(alpha, mu, gamma_1p1oa, z, n_bids){
  alpha_pv_1 = n_bids*alpha
  alpha_pv_2 = (gamma_1p1oa/mu)^alpha
  alpha_pv_3 = z^(alpha - 1)
  alpha_pv_4 = exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  deriv_alpha_pv_1 = n_bids
  deriv_alpha_pv_2 = (gamma_1p1oa/mu)^alpha*(log(gamma_1p1oa/mu) - digamma(1+1/alpha)/alpha)
  deriv_alpha_pv_3 = z^(alpha - 1)*log(z)
  deriv_alpha_pv_4 = -exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)*n_bids*(gamma_1p1oa/mu*z)^alpha*(
    log(gamma_1p1oa/mu*z) - digamma(1+1/alpha)/alpha
  )

  deriv_pv_alpha = deriv_alpha_pv_1*alpha_pv_2*alpha_pv_3*alpha_pv_4 +
    alpha_pv_1*deriv_alpha_pv_2*alpha_pv_3*alpha_pv_4 +
    alpha_pv_1*alpha_pv_2*deriv_alpha_pv_3*alpha_pv_4 +
    alpha_pv_1*alpha_pv_2*alpha_pv_3*deriv_alpha_pv_4

  return(deriv_pv_alpha)
}

auction__pv_alpha <- function(alpha, mu, z, n_bids){
  pv = n_bids*alpha*(gamma(1+1/alpha)/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma(1+1/alpha)/mu*z)^alpha)
}

grad(auction__pv_alpha,
     x = alpha,
     mu = mu,
     z = 2,
     n_bids = n_bids)


auction__deriv_pv_alpha(alpha = alpha,
                        mu = mu,
                        gamma_1p1oa = gamma_1p1oa,
                        z = 2,
                        n_bids = n_bids)



### Overall likelihood - checks out
f__w_integrand_z_fast <- function(z, w_bid, n_bids, mu, alpha, id_distrib, distrib_std_dev){

  b_z = f__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha)

  deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                     id_distrib = id_distrib,
                                     w_bid = w_bid,
                                     b_z = b_z)
  f_unobs = deriv_unobs[[1]]

  listFuncCall$argList$x = w_bid/b_z

  vals = n_bids*alpha*(gamma(1+1/alpha)/mu)^alpha*z^(alpha-1)*
    exp(-n_bids*(gamma(1+1/alpha)/mu*z)^alpha)*
    1/b_z*
    f_unobs

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}


auction__deriv_alpha_integrand <- function(alpha, w_bid, mu, gamma_1p1oa, z, n_bids, distrib_std_dev, id_distrib){
  b_z = vf__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)

  deriv_pv_alpha = auction__deriv_pv_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = n_bids)
  deriv_bid_alpha = auction__deriv_bid_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = n_bids)

  deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                     id_distrib = id_distrib,
                                     w_bid = w_bid,
                                     b_z = b_z)
  deriv_unobs_u = deriv_unobs[[2]]
  f_unobs = deriv_unobs[[1]]

  pv = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  vals = deriv_pv_alpha*1/b_z*f_unobs -
    pv*1/((b_z)^2)*deriv_bid_alpha*f_unobs -
    pv*1/b_z*deriv_unobs_u*w_bid/((b_z)^2)*deriv_bid_alpha

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0

  return(vals)
}

grad(f__w_integrand_z_fast,
     x = alpha,
     z = 2,
     n_bids = n_bids,
     mu = mu,
     w_bid = w_bid,
     id_distrib = 2,
     distrib_std_dev = 0.6)

auction__deriv_alpha_integrand(alpha = alpha,
                               w_bid = w_bid,
                               mu = mu,
                               gamma_1p1oa = gamma_1p1oa,
                               z = 2,
                               n_bids = n_bids,
                               distrib_std_dev = 0.6,
                               id_distrib = 2)


####MU
###Bid function - checks out

auction__deriv_bid_mu_integrand <- function(n_bids, gamma_1p1oa, alpha, mu, xi){
  val = exp(-(n_bids - 1)*(gamma_1p1oa/mu*xi)^alpha)*(n_bids - 1)*(gamma_1p1oa*xi)^alpha*
    alpha*mu^(-alpha - 1)
  return(val)
}

auction__deriv_bid_mu_integrate <- function(n_bids, gamma_1p1oa, alpha, mu, z){
  int = rep(NA, length(z))
  for (i in 1: length(z)){
    int[i]  = stats::integrate(auction__deriv_bid_mu_integrand,
                               n_bids = n_bids,
                               gamma_1p1oa = gamma_1p1oa,
                               alpha = alpha,
                               mu = mu,
                               lower = z[i],
                               upper = Inf,
                               abs.tol = 1e-10)$value
  }
  return(int)
}

auction__deriv_bid_mu <- function(z, alpha, mu, gamma_1p1oa, n_bids){
  mu_bid_1 = 1/alpha*mu/gamma_1p1oa*(n_bids - 1)^(-1/alpha)*pgamma((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha, 1/alpha, lower = FALSE)*gamma(1/alpha)
  mu_bid_2 = exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)

  deriv_mu_bid_1 = auction__deriv_bid_mu_integrate(n_bids = n_bids,
                                                   gamma_1p1oa = gamma_1p1oa,
                                                   alpha = alpha,
                                                   mu = mu,
                                                   z = z)


  deriv_mu_bid_2 = -exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)*(n_bids-1)*(gamma_1p1oa*z)^alpha*alpha*mu^(-alpha-1)

  deriv_bid_mu = deriv_mu_bid_1*mu_bid_2+
    mu_bid_1*deriv_mu_bid_2

  return(deriv_bid_mu)
}

f__bid_function_fast = function(cost, n_bids, mu, alpha){

  if (exp(-(n_bids-1)*(1/(mu/gamma(1+alpha))*cost)^alpha) == 0) {
    return(cost + mu/alpha*(n_bids-1)^(-1/alpha)*1/gamma(1+1/alpha)*
             ((n_bids-1)*(gamma(1+1/alpha)/mu*cost)^alpha)^(1/alpha-1))
  }

  cost + 1/alpha*(mu/gamma(1+1/alpha))*(n_bids-1)^(-1/alpha)*
    stats::pgamma((n_bids-1)*(1/(mu/gamma(1+1/alpha))*cost)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(n_bids-1)*(1/(mu/gamma(1+1/alpha))*cost)^alpha)
}

grad(f__bid_function_fast,
     x = mu,
     cost = 2,
     n_bids = n_bids,
     alpha = alpha)

auction__deriv_bid_mu(z = 2,
                      alpha = alpha,
                      mu = mu,
                      gamma_1p1oa = gamma_1p1oa,
                      n_bids = n_bids)


###Private values - checks out

auction__deriv_pv_mu <- function(alpha, mu, gamma_1p1oa, z, n_bids){
  mu_pv_1 = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)
  mu_pv_2 = exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  deriv_mu_pv_1 = -n_bids*alpha^2*(gamma_1p1oa)^alpha*mu^(-alpha - 1)*z^(alpha - 1)
  deriv_mu_pv_2 = exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)*alpha*n_bids*(gamma_1p1oa*z/mu)^alpha/mu

  deriv_pv_mu = deriv_mu_pv_1*mu_pv_2 + mu_pv_1*deriv_mu_pv_2
  return(deriv_pv_mu)
}

auction__pv_mu <- function(alpha, mu, gamma_1p1oa, z, n_bids){
  pv = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)
  return(pv)
}

grad(auction__pv_mu,
     x = mu,
     alpha = alpha,
     gamma_1p1oa = gamma_1p1oa,
     z = 2,
     n_bids = n_bids)

auction__deriv_pv_mu(alpha = alpha,
                     mu = mu,
                     gamma_1p1oa = gamma_1p1oa,
                     z = 2,
                     n_bids = n_bids)


### Overall likelihood - checks out

auction__deriv_mu_integrand <- function(alpha, w_bid, mu, gamma_1p1oa, z, n_bids, distrib_std_dev, id_distrib){
  b_z = f__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha)

  deriv_pv_mu = auction__deriv_pv_mu(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = n_bids)
  deriv_bid_mu = auction__deriv_bid_mu(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = n_bids)

  deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                     id_distrib = id_distrib,
                                     w_bid = w_bid,
                                     b_z = b_z)
  deriv_unobs_u = deriv_unobs[[2]]
  f_unobs = deriv_unobs[[1]]

  pv = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  vals = deriv_pv_mu*1/b_z*f_unobs -
    pv*1/((b_z)^2)*deriv_bid_mu*f_unobs -
    pv*1/b_z*deriv_unobs_u*w_bid/((b_z)^2)*deriv_bid_mu

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0

  return(vals)
}

f__w_integrand_z_fast <- function(z, w_bid, n_bids, mu, alpha, id_distrib, distrib_std_dev){

  b_z = f__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha)

  deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                     id_distrib = id_distrib,
                                     w_bid = w_bid,
                                     b_z = b_z)
  f_unobs = deriv_unobs[[1]]

  listFuncCall$argList$x = w_bid/b_z

  vals = n_bids*alpha*(gamma(1+1/alpha)/mu)^alpha*z^(alpha-1)*
    exp(-n_bids*(gamma(1+1/alpha)/mu*z)^alpha)*
    1/b_z*
    f_unobs

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}

grad(f__w_integrand_z_fast,
     x = mu,
     z = 2,
     n_bids = n_bids,
     alpha = alpha,
     w_bid = w_bid,
     id_distrib = 2,
     distrib_std_dev = 0.6)

auction__deriv_mu_integrand(alpha = alpha,
                               w_bid = w_bid,
                               mu = mu,
                               gamma_1p1oa = gamma_1p1oa,
                               z = 2,
                               n_bids = n_bids,
                               distrib_std_dev = 0.6,
                               id_distrib = 2)






