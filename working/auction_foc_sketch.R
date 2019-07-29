# Sketch implementation of analytical first derivatives for standard errors
##########
#UH

auction__deriv_unobs <- function(distrib_std_dev,
                                 id_distrib,
                                 w_bid,
                                 b_z){
  if (id_distrib == 1) {
    # dgamma
    deriv_unobs = auction__f_unobs_gamma(distrib_std_dev = distrib_std_dev,
                                         w_bid = w_bid,
                                         b_z = b_z)
  } else if (id_distrib == 2) {
    # dlnorm
    deriv_unobs = auction__f_unobs_lognorm(distrib_std_dev = distrib_std_dev,
                                           w_bid = w_bid,
                                           b_z = b_z)
  } else if (id_distrib == 3) {
    # dweibull
    deriv_unobs = auction__f_unobs_weibull(distrib_std_dev = distrib_std_dev,
                                           w_bid = w_bid,
                                           b_z = b_z)
  } else {
    res = list()
    res['err_code'] = 3
    res['err_msg'] = paste("Unknown id_distrib '", id_distrib, "'", sep='')
    auction__gen_err_msg(res)
  }
  return(deriv_unobs)
}

auction__f_unobs_gamma <- function(distrib_std_dev, w_bid, b_z){

  deriv_unobs_sig = -(1/distrib_std_dev^2)^(2+(1/distrib_std_dev^2))*(w_bid/b_z)^(1/distrib_std_dev^2-1)*
    exp(-(1/distrib_std_dev^2)*w_bid/b_z)*w_bid/b_z*(
      log(1/distrib_std_dev^2) + 1 - digamma(1/distrib_std_dev^2) + log(w_bid/b_z) - w_bid/b_z
    )

  deriv_unobs_u = -1/gamma(1/distrib_std_dev^2)*(1/distrib_std_dev^2)^(1+1/distrib_std_dev^2)*
    (w_bid/b_z)^(1/distrib_std_dev^2 -2)*(distrib_std_dev^2-1+w_bid/b_z)*
    exp(-1/distrib_std_dev^2 * w_bid/b_z)

  f_unobs = (1/(distrib_std_dev^2))^(1/(distrib_std_dev^2))/gamma(1/(distrib_std_dev^2))*
    (w_bid/b_z)^(1/(distrib_std_dev^2)-1)*exp(-1/(distrib_std_dev^2)*w_bid/b_z)
  return(list(f_unobs, deriv_unobs_u, deriv_unobs_sig))
}

auction__f_unobs_lognorm <- function(distrib_std_dev, w_bid, b_z){

  deriv_unobs_sig = -1/(w_bid/b_z*sqrt(2*pi))*1/(2*(log(1+distrib_std_dev^2))^(3/2))*1/(1+distrib_std_dev^2) *
    exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))
  -1/(w_bid/b_z*sqrt(2*pi*log(1+distrib_std_dev^2)))*exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))*
  (
    2*(log(1+distrib_std_dev^2)-1)*(log(w_bid/b_z)-1/2*log(1+distrib_std_dev^2))/(4*(log(1+distrib_std_dev^2))^2*(1+distrib_std_dev^2))
  )

  deriv_unobs_u = exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))*
    (2*log(w_bid/b_z)+3*log(1+distrib_std_dev^2))/(2*(w_bid/b_z)^2*(log(1+distrib_std_dev^2))^(3/2)*sqrt(2*pi))

  f_unobs = 1/(w_bid/b_z*sqrt(2*pi*log(1+distrib_std_dev^2)))*
    exp(-(log(w_bid/b_z)+1/2*log(1+distrib_std_dev^2))^2/(2*log(1+distrib_std_dev^2)))
  return(list(f_unobs, deriv_unobs_u, deriv_unobs_sig))
}

auction__f_unobs_weibull <- function(distrib_std_dev, w_bid, b_z){

  listParam = auction__get_unobs_params(distrib_std_dev, 3)
  shape = listParam$shape

  dshape = (gamma(1+1/shape))^2*(
    4/(shape^4)*gamma(1+2/shape)*digamma(1+2/shape)*gamma(1+1/shape) - 1/(shape^4)*gamma(1+2/shape)*gamma(1+1/shape)*digamma(1+1/shape)
  )

  df = gamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  - 1/shape*gamma(1+1/shape)*digamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  + shape*gamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(shape - 1)*(
    log(w_bid/b_z*gamma(1+1/shape)) - (shape - 1)*digamma(1+1/shape)/(shape^2)
  )*exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  - 1/shape*gamma(1+1/shape)*((w_bid/b_z)*gamma(1+1/shape))^(2*shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)*(
      log(w_bid/b_z*gamma(1+1/shape)) - digamma(1+1/shape)/shape
    )

  deriv_unobs_sig = df*dshape

  deriv_unobs_u = shape*(gamma(1+1/shape))^shape*exp(-(w_bid/b_z*gamma(1+1/shape))^shape)*(
    (shape - 1)*(w_bid/b_z)^(shape -2) - (w_bid/b_z)^shape*shape*(gamma(1+1/shape))^shape
  )

  f_unobs = shape*gamma(1+1/shape)*(w_bid/b_z*gamma(1+1/shape))^(shape - 1)*
    exp(-(w_bid/b_z*gamma(1+1/shape))^shape)
  return(list(f_unobs, deriv_unobs_u, deriv_unobs_sig))
}

############

#alpha
auction__deriv_bid_alpha_integrand <- function(n_bids, gamma_1p1oa, alpha, mu, xi){
  val = -exp(-(n_bids - 1)*(gamma_1p1oa/mu*xi)^alpha)*(n_bids - 1)*(gamma_1p1oa/mu*xi)^alpha*(
    log(gamma_1p1oa/mu*xi) - digamma(1+1/alpha)/alpha
  )
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
  alpha_bid_1 = 1/alpha
  alpha_bid_2 = mu/gamma_1p1oa
  alpha_bid_3 = (n_bids - 1)^(-1/alpha)
  alpha_bid_4 = pgamma((n_bids- 1)*(gamma_1p1oa/mu*z)^alpha, 1/alpha, lower = FALSE)*gamma(1/alpha)
  alpha_bid_5 = exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)


  deriv_alpha_bid_1 = -1/alpha^2
  deriv_alpha_bid_2 = mu/(alpha^2) *digamma(1+1/alpha)/(gamma_1p1oa)
  deriv_alpha_bid_3 = (n_bids - 1)^(-1/alpha)*log(n_bids - 1)/(alpha^2)
  deriv_alpha_bid_4 = auction__deriv_bid_alpha_integrate(n_bids = n_bids,
                                                         gamma_1p1oa = gamma_1p1oa,
                                                         alpha = alpha,
                                                         mu = mu,
                                                         z = z)
  deriv_alpha_bid_5 = exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)*
    (n_bids - 1)*(gamma_1p1oa/mu*z)^alpha*(
      log(gamma_1p1oa/mu*z) - digamma(1+1/alpha)/alpha
    )

  deriv_bid_alpha = deriv_alpha_bid_1*alpha_bid_2*alpha_bid_3*alpha_bid_4*alpha_bid_5 +
    alpha_bid_1*deriv_alpha_bid_2*alpha_bid_3*alpha_bid_4*alpha_bid_5 +
    alpha_bid_1*alpha_bid_2*deriv_alpha_bid_3*alpha_bid_4*alpha_bid_5 +
    alpha_bid_1*alpha_bid_2*alpha_bid_3*deriv_alpha_bid_4*alpha_bid_5 +
    alpha_bid_1*alpha_bid_2*alpha_bid_3*alpha_bid_4*deriv_alpha_bid_5

  return(deriv_bid_alpha)
}

auction__vec_deriv_bid_alpha = Vectorize(auction__deriv_bid_alpha, vectorize.args = "z")

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


auction__deriv_alpha_integrate <- function(data_vec){
  val = stats::integrate(auction__deriv_alpha_integrand,
                         w_bid=data_vec[1],
                         n_bids=data_vec[2],
                         mu=data_vec[3],
                         alpha=data_vec[4],
                         gamma_1p1oa=data_vec[5],
                         distrib_std_dev=data_vec[6],
                         id_distrib=data_vec[7],
                         lower = 0,
                         upper = Inf,
                         abs.tol = 1e-10)

  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

auction__deriv_ll_alpha <- function(x0, dat_X, dat__winning_bid, dat__n_bids, cl, listFuncCall){
  listIdx = auction__x0_indices()

  if(x0[listIdx$unobs_dist_param] <= 0.1) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_mu] <= 0) > 0) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_a] <= 0.01) > 0) return(-Inf)

  param.h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1)]
  v.h = exp( colSums(param.h * dat_X) )

  listFuncCall$argList = auction__get_unobs_params(
    distrib_std_dev = x0[listIdx$unobs_dist_param],
    id_distrib = listFuncCall$funcID)

  v.f_w = parallel::parApply(cl = cl,
                             X = cbind(dat__winning_bid/v.h,
                                       dat__n_bids,
                                       x0[listIdx$pv_weibull_mu],
                                       x0[listIdx$pv_weibull_a],
                                       gamma(1 + 1/x0[listIdx$pv_weibull_a]) ),
                             MARGIN = 1,
                             FUN = f__funk,
                             listFuncCall=listFuncCall)

  v.deriv_f_w_alpha = parallel::parApply(cl = cl,
                                         X = cbind(dat__winning_bid/v.h,
                                                   dat__n_bids,
                                                   x0[listIdx$pv_weibull_mu],
                                                   x0[listIdx$pv_weibull_a],
                                                   gamma(1 + 1/x0[listIdx$pv_weibull_a]),
                                                   x0[listIdx$unobs_dist_param],
                                                   listFuncCall$funcID),
                                         MARGIN = 1,
                                         FUN = auction__deriv_alpha_integrate)
  v.f_deriv_y_alpha = 1/(v.f_w*v.h)*v.deriv_f_w_alpha
  deriv_alpha_log_likelihood = -sum(v.f_deriv_y_alpha)
  return(list(deriv_alpha_log_likelihood,v.f_deriv_y_alpha))
}

##########
## mu
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
  mu_bid_1 = 1/alpha*mu/gamma_1p1oa*(n_bids - 1)^(-1/alpha)
  mu_bid_2 = pgamma((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha, 1/alpha, lower = FALSE)*gamma(1/alpha)
  mu_bid_3 = exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)

  deriv_mu_bid_1 = 1/(alpha*gamma_1p1oa)*(n_bids-1)^(-1/alpha)
  deriv_mu_bid_2 = auction__deriv_bid_mu_integrate(n_bids = n_bids,
                                                   gamma_1p1oa = gamma_1p1oa,
                                                   alpha = alpha,
                                                   mu = mu,
                                                   z = z)


  deriv_mu_bid_3 = -exp((n_bids - 1)*(gamma_1p1oa/mu*z)^alpha)*(n_bids-1)*(gamma_1p1oa*z)^alpha*alpha*mu^(-alpha-1)

  deriv_bid_mu = deriv_mu_bid_1*mu_bid_2*mu_bid_3+
    mu_bid_1*deriv_mu_bid_2*mu_bid_3+
    mu_bid_1*mu_bid_2*deriv_mu_bid_3

  return(deriv_bid_mu)
}

auction__vec_deriv_bid_mu = Vectorize(auction__deriv_bid_mu, vectorize.args = "z")

auction__deriv_pv_mu <- function(alpha, mu, gamma_1p1oa, z, n_bids){
  mu_pv_1 = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)
  mu_pv_2 = exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  deriv_mu_pv_1 = -n_bids*alpha^2*(gamma_1p1oa)^alpha*mu^(-alpha - 1)*z^(alpha - 1)
  deriv_mu_pv_2 = exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)*alpha*n_bids*(gamma_1p1oa*z)^alpha*mu*(-alpha-1)

  deriv_pv_mu = deriv_mu_pv_1*mu_pv_2 + mu_pv_1*deriv_mu_pv_2
  return(deriv_pv_mu)
}

auction__deriv_mu_integrand <- function(alpha, w_bid, mu, gamma_1p1oa, z, n_bids, distrib_std_dev, id_distrib){
  b_z = vf__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)

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

auction__deriv_mu_integrate <- function(data_vec){
  val = stats::integrate(auction__deriv_mu_integrand,
                         w_bid=data_vec[1],
                         n_bids=data_vec[2],
                         mu=data_vec[3],
                         alpha=data_vec[4],
                         gamma_1p1oa=data_vec[5],
                         distrib_std_dev=data_vec[6],
                         id_distrib=data_vec[7],
                         lower = 0,
                         upper = Inf,
                         abs.tol = 1e-10)
  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

auction__deriv_ll_mu <- function(x0, dat__winning_bid, dat__n_bids, dat_X, cl, listFuncCall){
  listIdx = auction__x0_indices()

  if(x0[listIdx$unobs_dist_param] <= 0.1) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_mu] <= 0) > 0) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_a] <= 0.01) > 0) return(-Inf)

  param.h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1)]
  v.h = exp( colSums(param.h * dat_X) )

  listFuncCall$argList = auction__get_unobs_params(
    distrib_std_dev = x0[listIdx$unobs_dist_param],
    id_distrib = listFuncCall$funcID)

  v.f_w = parallel::parApply(cl = cl,
                             X = cbind(dat__winning_bid/v.h,
                                       dat__n_bids,
                                       x0[listIdx$pv_weibull_mu],
                                       x0[listIdx$pv_weibull_a],
                                       gamma(1 + 1/x0[listIdx$pv_weibull_a]) ),
                             MARGIN = 1,
                             FUN = f__funk,
                             listFuncCall=listFuncCall)

  v.deriv_f_w_mu = parallel::parApply(cl = cl,
                                      X = cbind(dat__winning_bid/v.h,
                                                dat__n_bids,
                                                x0[listIdx$pv_weibull_mu],
                                                x0[listIdx$pv_weibull_a],
                                                gamma(1 + 1/x0[listIdx$pv_weibull_a]),
                                                x0[listIdx$unobs_dist_param],
                                                listFuncCall$funcID),
                                      MARGIN = 1,
                                      FUN = auction__deriv_mu_integrate)
  v.f_deriv_y_mu = 1/(v.f_w*v.h)*v.deriv_f_w_mu
  deriv_mu_log_likelihood = -sum(v.f_deriv_y_mu)
  return(list(deriv_mu_log_likelihood,v.f_deriv_y_mu))
}

######
#sigma^2
auction__deriv_sig_integrand <- function(distrib_std_dev, w_bid, id_distrib, z, n_bids, mu, alpha, gamma_1p1oa){
  b_z = vf__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)
  deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                     id_distrib = id_distrib,
                                     w_bid = w_bid,
                                     b_z = b_z)
  deriv_unobs_sig = deriv_unobs[[3]]
  f_unobs = deriv_unobs[[1]]

  pv = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  vals = pv*1/b_z*deriv_unobs_sig

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0

  return(vals)

}

auction__deriv_sig_integrate <- function(data_vec){
  val = stats::integrate(auction__deriv_sig_integrand,
                         w_bid=data_vec[1],
                         n_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                         gamma_1p1oa=data_vec[5],
                         distrib_std_dev=data_vec[6],
                         id_distrib=data_vec[7],
                         lower = 0,
                         upper = Inf,
                         abs.tol = 1e-10)
  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

auction__deriv_ll_sig <- function(x0, dat__winning_bid, dat__n_bids, dat_X, cl, listFuncCall){
  listIdx = auction__x0_indices()

  if(x0[listIdx$unobs_dist_param] <= 0.1) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_mu] <= 0) > 0) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_a] <= 0.01) > 0) return(-Inf)

  param.h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1)]
  v.h = exp( colSums(param.h * dat_X) )

  listFuncCall$argList = auction__get_unobs_params(
    distrib_std_dev = x0[listIdx$unobs_dist_param],
    id_distrib = listFuncCall$funcID)

  v.f_w = parallel::parApply(cl = cl,
                             X = cbind(dat__winning_bid/v.h,
                                       dat__n_bids,
                                       x0[listIdx$pv_weibull_mu],
                                       x0[listIdx$pv_weibull_a],
                                       gamma(1 + 1/x0[listIdx$pv_weibull_a]) ),
                             MARGIN = 1,
                             FUN = f__funk,
                             listFuncCall=listFuncCall)
  v.deriv_f_w_sig = parallel::parApply(cl = cl,
                                       X = cbind(dat__winning_bid/v.h,
                                                 dat__n_bids,
                                                 x0[listIdx$pv_weibull_mu],
                                                 x0[listIdx$pv_weibull_a],
                                                 gamma(1 + 1/x0[listIdx$pv_weibull_a]),
                                                 x0[listIdx$unobs_dist_param],
                                                 listFuncCall$funcID),
                                       MARGIN = 1,
                                       FUN = auction__deriv_sig_integrate)
  v.f_deriv_y_sig = 1/(v.f_w*v.h)*v.deriv_f_w_sig
  deriv_sig_log_likelihood = -sum(v.f_deriv_y_sig)
  return(list(deriv_sig_log_likelihood, v.f_deriv_y_sig))
}
#######
#betas
auction__deriv_beta_integrand <- function(w_bid, mu, alpha, gamma_1p1oa, distrib_std_dev, id_distrib, n_bids, x_k, z){

  b_z = auctionmodel:::vf__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)
  deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                     id_distrib = id_distrib,
                                     w_bid = w_bid,
                                     b_z = b_z)
  deriv_unobs_u = deriv_unobs[[2]]
  f_unobs = deriv_unobs[[1]]

  pv = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)

  deriv_h = -pv*1/((b_z)^2)*deriv_unobs_u*w_bid*x_k

  deriv_h[(gamma_1p1oa/mu)^alpha == Inf] = 0
  deriv_h[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0

  vals = deriv_h
  return(vals)
}

auction__deriv_beta_integrate <- function(data_vec){
  val = stats::integrate(auction__deriv_beta_integrand,
                         w_bid=data_vec[1],
                         n_bids=data_vec[2],
                         mu=data_vec[3],
                         alpha=data_vec[4],
                         gamma_1p1oa=data_vec[5],
                         distrib_std_dev=data_vec[6],
                         id_distrib=data_vec[7],
                         x_k = data_vec[8],
                         lower = 0,
                         upper = Inf,
                         abs.tol = 1e-10)
  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

auction__deriv_ll_beta <- function(x0, dat__winning_bid, dat__n_bids, dat_X, cl, listFuncCall, x_k){
  listIdx = auctionmodel:::auction__x0_indices()

  if(x0[listIdx$unobs_dist_param] <= 0.1) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_mu] <= 0) > 0) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_a] <= 0.01) > 0) return(-Inf)

  param.h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1)]
  v.h = exp( colSums(param.h * dat_X) )

  listFuncCall$argList = auctionmodel:::auction__get_unobs_params(
    distrib_std_dev = x0[listIdx$unobs_dist_param],
    id_distrib = listFuncCall$funcID)

  v.f_w = parallel::parApply(cl = cl,
                             X = cbind(dat__winning_bid/v.h,
                                       dat__n_bids,
                                       x0[listIdx$pv_weibull_mu],
                                       x0[listIdx$pv_weibull_a],
                                       gamma(1 + 1/x0[listIdx$pv_weibull_a]) ),
                             MARGIN = 1,
                             FUN = auctionmodel:::f__funk,
                             listFuncCall=listFuncCall)

  v.deriv_f_w_beta = parallel::parApply(cl = cl,
                                        X = cbind(dat__winning_bid/v.h,
                                                  dat__n_bids,
                                                  x0[listIdx$pv_weibull_mu],
                                                  x0[listIdx$pv_weibull_a],
                                                  gamma(1 + 1/x0[listIdx$pv_weibull_a]),
                                                  x0[listIdx$unobs_dist_param],
                                                  listFuncCall$funcID,
                                                  x_k),
                                        MARGIN = 1,
                                        FUN = auction__deriv_beta_integrate
  )

  v.f_deriv_y_beta = 1/(v.f_w*v.h)*v.deriv_f_w_beta
  deriv_beta_log_likelihood = -sum(v.f_deriv_y_beta)
  auction__deriv_ll_beta = list(deriv_beta_log_likelihood, v.f_deriv_y_beta)
  return(auction__deriv_ll_beta)
}

auction__deriv_ll_beta_wrap <- function(x0, dat__winning_bid, dat__n_bids, dat_X, cl, listFuncCall){
  deriv_betas = list()

  beta_length = length(x0)-listIdx$x_terms__start + 1
  for (k in 1:beta_length){
    nam = paste("x_",k,sep="")
    assign(nam,dat_X[k,])
    deriv_betas[[k]] = auction__deriv_ll_beta(x0 = x0,
                                              dat__winning_bid = dat__winning_bid,
                                              dat__n_bids = dat__n_bids,
                                              dat_X = dat_X,
                                              cl = cl,
                                              listFuncCall = listFuncCall,
                                              x_k = x_k)
  }
  return(deriv_betas)
}


#######
#standard errors
auction__se <- function(x0, dat__winning_bid, dat__n_bids, dat_X, cl, listFuncCall, obs){
  auction__deriv_ll_alpha = auction__deriv_ll_alpha(x0 = x0,
                                                    dat__winning_bid = dat__winning_bid,
                                                    dat__n_bids = dat__n_bids,
                                                    dat_X = dat_X,
                                                    cl = cl,
                                                    listFuncCall = listFuncCall)
  v.f_deriv_y_alpha = auction__deriv_ll_alpha[[2]]

  auction__deriv_ll_mu = auction__deriv_ll_mu(x0 = x0,
                                              dat__winning_bid = dat__winning_bid,
                                              dat__n_bids = dat__n_bids,
                                              dat_X = dat_X,
                                              cl = cl,
                                              listFuncCall = listFuncCall)
  v.f_deriv_y_mu = auction__deriv_ll_mu[[2]]

  auction__deriv_ll_beta = auction__deriv_ll_beta_wrap(x0 = x0,
                                                       dat__winning_bid = dat__winning_bid,
                                                       dat__n_bids = dat__n_bids,
                                                       dat_X = dat_X,
                                                       cl = cl,
                                                       listFuncCall = listFuncCall)


  auction__deriv_ll_sig = auction__deriv_ll_sig(x0 = x0,
                                                dat__winning_bid = dat__winning_bid,
                                                dat__n_bids = dat__n_bids,
                                                dat_X = dat_X,
                                                cl = cl,
                                                listFuncCall = listFuncCall)
  v.f_deriv_y_sig = auction__deriv_ll_sig[[2]]

  auction__deriv_ll_beta = auction__deriv_ll_beta_wrap(x0 = x0,
                                                       dat__winning_bid = dat__winning_bid,
                                                       dat__n_bids = dat__n_bids,
                                                       dat_X = dat_X,
                                                       cl = cl,
                                                       listFuncCall = listFuncCall)


  beta_length = length(x0)-listIdx$x_terms__start + 1
  for (k in 1:beta_length){
    nam = paste("v.f_deriv_y_beta_",k,sep="")
    assign(nam,auction__deriv_ll_beta[[k]][[2]])
  }

  all_prod = list()
  score = list()
  for (i in 1:obs){
    score_i = c(v.f_deriv_y_alpha[i], v.f_deriv_y_mu[i], v.f_deriv_y_sig[i])
    for (k in 1:beta_length){
      score_i[3+k] = v.f_deriv_y_beta_k[i]
    }
    prod_i = score_i %*% t(score_i)
    all_prod[[i]] = prod_i
    score[[i]] = score_i
  }

  sum_prodmat = Reduce('+',all_prod)

  info_mat = 1/obs*sum_prodmat

  varcov_params = solve(info_mat)

  var_params = diag(varcov_params)

  se_params = sqrt(var_params)

  return(se_params)
}

