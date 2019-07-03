rm(list = ls())
library(auctionmodel)
library(devtools)
library(testthat)

context("Bid function, log likelihood")

test_that("Bid function",{
  n_bids = sample(2:10, 1, replace = TRUE)
  mu = 2
  alpha = 3
  gamma_1p1oa = gamma(1 + 1/alpha)
  cost = rep(NA, n_bids)
  # Generate cost
  cost = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(n_bids)))^(1/alpha)
  bid = rep(NA, n_bids)

  for (i in 1:n_bids){
  bid[i] = auctionmodel:::f__bid_function_fast(cost = cost[i],
                                             n_bids = n_bids,
                                             mu = mu,
                                             alpha = alpha,
                                             gamma_1p1oa = gamma_1p1oa)

  expect_equal(bid[i],
               cost[i] + 1/alpha*(mu/gamma_1p1oa)*(n_bids-1)^(-1/alpha)*
               stats::pgamma((n_bids-1)*(1/(mu/gamma_1p1oa)*cost[i])^alpha, 1/alpha, lower=FALSE)*
               gamma(1/alpha)*
               1/exp(-(n_bids-1)*(1/(mu/gamma_1p1oa)*cost[i])^alpha))
  }
  expect_gte(min(bid),0)
  expect_equal(is.numeric(bid),
               TRUE)

  # test integrand at zero
  alpha = 0
  gamma_1p1oa = gamma(1 + 1/alpha)
  cost = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(n_bids)))^(1/alpha)
  for (i in 1:n_bids){
    bid[i] = auctionmodel:::f__bid_function_fast(cost = cost[i],
                                                 n_bids = n_bids,
                                                 mu = mu,
                                                 alpha = alpha,
                                                 gamma_1p1oa = gamma_1p1oa)

    expect_equal(bid[i],
                 cost[i] + mu/alpha*(n_bids-1)^(-1/alpha)*1/gamma_1p1oa*
                 ((n_bids-1)*(gamma_1p1oa/mu*cost[i])^alpha)^(1/alpha-1))
  }


  # test vectorize
  alpha = 3
  gamma_1p1oa = gamma(1 + 1/alpha)
  n_bids = sample(2:10, 1, replace = TRUE)
  mu = 5
  cost = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(n_bids)))^(1/alpha)
  bid = rep(NA, n_bids)
  for (i in 1:n_bids){
    bid[i] = auctionmodel:::f__bid_function_fast(cost = cost[i],
                                                 n_bids = n_bids,
                                                 mu = mu,
                                                 alpha = alpha,
                                                 gamma_1p1oa = gamma_1p1oa)
  }

  bid1 = auctionmodel:::vf__bid_function_fast(cost = cost,
                                              n_bids = n_bids,
                                              mu = mu,
                                              alpha = alpha,
                                              gamma_1p1oa = gamma_1p1oa)
  expect_equal(bid, bid1)
})


test_that("Observed bid density integrand", {
  n_bids = sample(2:10, 1, replace = TRUE)
  mu = 2
  alpha = 3
  gamma_1p1oa = gamma(1 + 1/alpha)
  cost = rep(NA, n_bids)
  # Generate cost
  cost = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(n_bids)))^(1/alpha)
  # Get winning bid
  w_cost = min(cost)
  w_bid = auctionmodel:::f__bid_function_fast(cost = w_cost,
                                              n_bids = n_bids,
                                              mu = mu,
                                              alpha = alpha,
                                              gamma_1p1oa = gamma_1p1oa)
  # Get u dist + parameters
  u_dist = "dgamma"
  for (funcName in u_dist){
  sFuncName = as.character(funcName)
  }

  listFuncCall = list(funcName = sFuncName,
                      funcID = auctionmodel:::auction__get_id_distrib(sFuncName = sFuncName))

  sigma = .9
  listFuncCall$argList = auctionmodel:::auction__get_unobs_params(
    distrib_std_dev = sigma,
    id_distrib = listFuncCall$funcID)

  # integrand
  vf.w_int = auctionmodel:::vf__w_integrand_z_fast(z = cost,
                                        w_bid = w_bid,
                                        n_bids = n_bids,
                                        mu = mu,
                                        alpha = alpha,
                                        gamma_1p1oa = gamma_1p1oa,
                                        listFuncCall = listFuncCall)

  # test generic return
  b_z = auctionmodel:::vf__bid_function_fast(cost=cost, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  listFuncCall$argList$x = w_bid/b_z

  expect_equal(vf.w_int, n_bids*alpha*(gamma(1 + 1/alpha)/mu)^alpha*cost^(alpha-1)*
                 exp(-n_bids*(gamma(1 + 1/alpha)/mu*cost)^alpha)*
                 1/b_z*
                 do.call(
                   match.fun(listFuncCall$funcName), # UH term
                   listFuncCall$argList
                 ))

  # # limit cases ?
  # mu = 5
  # alpha = ?
  # gamma_1p1oa = gamma(1 + 1/alpha)
  # cost = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(n_bids)))^(1/alpha)
  # # Get winning bid
  # w_cost = min(cost)
  # w_bid = auctionmodel:::f__bid_function_fast(cost = w_cost,
  #                                             n_bids = n_bids,
  #                                             mu = mu,
  #                                             alpha = alpha,
  #                                             gamma_1p1oa = gamma_1p1oa)
  #
  # vf.w_int = auctionmodel:::vf__w_integrand_z_fast(z = cost,
  #                                                  w_bid = w_bid,
  #                                                  n_bids = n_bids,
  #                                                  mu = mu,
  #                                                  alpha = alpha,
  #                                                  gamma_1p1oa = gamma_1p1oa,
  #                                                  listFuncCall = listFuncCall)
  #
  # expect_equal(vf.w_int, rep(0, n_bids))

  test_that("Observed bid denasity",{
    #setup
    n_bids = sample(2:10, 1, replace = TRUE)
    mu = 2
    alpha = 3
    gamma_1p1oa = gamma(1 + 1/alpha)
    cost = rep(NA, n_bids)
    cost = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(n_bids)))^(1/alpha)
    w_cost = min(cost)
    w_bid = auctionmodel:::f__bid_function_fast(cost = w_cost,
                                                n_bids = n_bids,
                                                mu = mu,
                                                alpha = alpha,
                                                gamma_1p1oa = gamma_1p1oa)
    u_dist = "dgamma"
    for (funcName in u_dist){
      sFuncName = as.character(funcName)
    }

    listFuncCall = list(funcName = sFuncName,
                        funcID = auctionmodel:::auction__get_id_distrib(sFuncName = sFuncName))

    sigma = .9
    listFuncCall$argList = auctionmodel:::auction__get_unobs_params(
      distrib_std_dev = sigma,
      id_distrib = listFuncCall$funcID)

    data_vec = data.frame(cbind(w_bid,
                                n_bids,
                                mu,
                                alpha,
                                gamma_1p1oa))

    # test
    v.f_w = auctionmodel:::f__funk(data_vec = data_vec,
                                   listFuncCall = listFuncCall)
  })



})

