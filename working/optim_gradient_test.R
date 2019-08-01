obs = 300
mu = 4
alpha = 3
sigma = 0.6

data = auctionmodel:::auction_generate_data(obs = obs,
                              mu = mu,
                              alpha = alpha,
                              sigma = sigma,
                              u_dist = "dgamma")



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

hTracker = auction__tracker__build(hEnv_tmp=hEnv_tmp, report=0)

stats::optim(par=x0,
             fn=f__ll_parallel,
             gr=auction__deriv_gradient,
             control=list(maxit = 2000, parscale = c(1, 0.1, 1, 0.1, 0.1)),
             dat__winning_bid=data$winning_bid,
             dat__n_bids=data$n_bids,
             dat_X=t( data[ ! names(data) %in% c("winning_bid", "n_bids") ] ),
             listFuncCall=listFuncCall,
             hTracker=hTracker,
             method = "BFGS",
             cl=cl)

stats::optim(par=x0,
             fn=f__ll_parallel,
             control=list(maxit = 2000, parscale = c(1, 0.1, 1, 0.1, 0.1)),
             dat__winning_bid=data$winning_bid,
             dat__n_bids=data$n_bids,
             dat_X=t( data[ ! names(data) %in% c("winning_bid", "n_bids") ] ),
             listFuncCall=listFuncCall,
             hTracker=hTracker,
             cl=cl)



