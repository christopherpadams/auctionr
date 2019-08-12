
#### test
x0 =  c(8, 2, .5, .4, 0)
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
                                  "auction__deriv_unobs",
                                  "auction__f_unobs_gamma",
                                  "auction__f_unobs_lognorm",
                                  "auction__f_unobs_weibull"),
                        envir = environment(auction_model))


x0 = c(7.93369419885821, 0.66080453704268, 0.748151129180016,
       0.399372886551224, -0.00062711344877579)
auction__deriv_ll_alpha(x0 = x0,
                      dat__winning_bid = dat__winning_bid,
                      dat__n_bids = dat__n_bids,
                      dat_X = dat_X,
                      cl = cl,
                      listFuncCall = listFuncCall)

data_vec = list()
for (i in 1:obs){
data_vec[[i]] = cbind(dat__winning_bid[i], dat__n_bids[i], 7.93369419885821, 3, gamma(1+1/3),  0.748151129180016, 1)
}

integrate = rep(NA, obs)
for (i in 1:300){
integrate[i] = auction__deriv_alpha_integrate(data_vec = data_vec[[i]])
}

for (i in 1:300){
  integrate[i] = auction__deriv_alpha_integrand(alpha = alpha,
                                                w_bid = dat__winning_bid[i],
                                                mu = mu,
                                                gamma_1p1oa = gamma(1+1/alpha),
                                                z = z,
                                                n_bids = 5,
                                                distrib_std_dev = distrib_std_dev,
                                                id_distrib = id_distrib)
}

auction__deriv_alpha_integrate <- function(data_vec){
  out = tryCatch({
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
    return(val$value)
  }, error = function(cond){
    val_e = tryCatch({
    val_1 = stats::integrate(auction__deriv_alpha_integrand,
                           w_bid=data_vec[1],
                           n_bids=data_vec[2],
                           mu=data_vec[3],
                           alpha=data_vec[4],
                           gamma_1p1oa=data_vec[5],
                           distrib_std_dev=data_vec[6],
                           id_distrib=data_vec[7],
                           lower = 0,
                           upper = 1,
                           abs.tol = 1e-10)
    val_2 = stats::integrate(auction__deriv_alpha_integrand,
                           w_bid=data_vec[1],
                           n_bids=data_vec[2],
                           mu=data_vec[3],
                           alpha=data_vec[4],
                           gamma_1p1oa=data_vec[5],
                           distrib_std_dev=data_vec[6],
                           id_distrib=data_vec[7],
                           lower = 1,
                           upper = Inf,
                           abs.tol = 1e-10)
    val_e = val_1$value + val_2$value
    return(val_e)
    }, error = function(cond){
      val_1 = stats::integrate(auction__deriv_alpha_integrand,
                               w_bid=data_vec[1],
                               n_bids=data_vec[2],
                               mu=data_vec[3],
                               alpha=data_vec[4],
                               gamma_1p1oa=data_vec[5],
                               distrib_std_dev=data_vec[6],
                               id_distrib=data_vec[7],
                               lower = 0,
                               upper = 1e-5,
                               abs.tol = 1e-10)
      val_2 = stats::integrate(auction__deriv_alpha_integrand,
                               w_bid=data_vec[1],
                               n_bids=data_vec[2],
                               mu=data_vec[3],
                               alpha=data_vec[4],
                               gamma_1p1oa=data_vec[5],
                               distrib_std_dev=data_vec[6],
                               id_distrib=data_vec[7],
                               lower = 1e-5,
                               upper = 1,
                               abs.tol = 1e-10)
      val_3 = stats::integrate(auction__deriv_alpha_integrand,
                               w_bid=data_vec[1],
                               n_bids=data_vec[2],
                               mu=data_vec[3],
                               alpha=data_vec[4],
                               gamma_1p1oa=data_vec[5],
                               distrib_std_dev=data_vec[6],
                               id_distrib=data_vec[7],
                               lower = 1,
                               upper = Inf,
                               abs.tol = 1e-10)
      val_e2 = val_1$value + val_2$value + val_3$value
      return(val_e2)
    })
    return(val_e)
  })
  return(out)
}






z = 1:100
plot(z, auction__deriv_alpha_integrand(alpha = alpha,
                               w_bid = 3.059029,
                               mu = mu,
                               gamma_1p1oa = gamma(1+1/alpha),
                               z = z,
                               n_bids = 5,
                               distrib_std_dev = distrib_std_dev,
                               id_distrib = id_distrib), ylab = "")



b_z = vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)

deriv_pv_alpha = auction__deriv_pv_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =cost, n_bids = 5)
deriv_bid_alpha = auction__deriv_bid_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =cost, n_bids = 5)

deriv_unobs = auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                   id_distrib = id_distrib,
                                   w_bid = w_bid,
                                   b_z = b_z)


deriv_unobs_u = deriv_unobs[[2]]
f_unobs = deriv_unobs[[1]]

pv = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha - 1)*exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)


z = 1:1000/100
plot(z,
     auction__deriv_pv_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z = z, n_bids = 5),
     ylab = "")


plot(z, auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                             id_distrib = id_distrib,
                             w_bid = w_bid,
                             b_z = vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))[[2]],
     ylab = "")


plot(z, auction__deriv_pv_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z = z, n_bids = 5)*
       1/vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa)*
       auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                            id_distrib = id_distrib,
                            w_bid = w_bid,
                            b_z = vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))[[1]],
     ylab="")


plot(z, 5*0.66*(gamma(1+1/0.66)/10)^0.66*z^(0.66 - 1)*exp(-5*(gamma(1+1/0.66)/10*z)^0.66)
*1/((vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))^2)*
  auction__deriv_bid_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = 5)*
  auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                                         id_distrib = id_distrib,
                                         w_bid = w_bid,
                                         b_z = vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))[[1]],
ylab = "")


plot(z, 5*0.66*(gamma(1+1/0.66)/10)^0.66*z^(0.66 - 1)*exp(-5*(gamma(1+1/0.66)/10*z)^0.66)
     *1/1/(vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))
     *auction__deriv_unobs(distrib_std_dev = distrib_std_dev,
                           id_distrib = id_distrib,
                           w_bid = w_bid,
                           b_z = vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))[[2]]
     *1.98/((vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa))^2)*
  auction__deriv_bid_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = 5),
  ylab = "")

# plot(z, auction__deriv_bid_alpha(alpha = alpha, mu = mu, gamma_1p1oa = gamma_1p1oa, z =z, n_bids = 5),
#      ylab = "")


# plot(z, vf__bid_function_fast(cost=z, n_bids=5, mu=mu, alpha=alpha, gamma_1p1oa = gamma_1p1oa),
#      ylab = "")



