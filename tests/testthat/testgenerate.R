rm(list = ls())
library(auctionmodel)
library(devtools)
library(testthat)

context("Generate data")

# tests without controls
test_that("Requires appropriate parameters, obs, nbids", {
  obs = 100
  mu = 5
  alpha = 3
  sigma = .6
  # Requires parameters
  expect_equal(auction_generate_data(obs = obs),
               NULL)
  # Requires obs
  expect_equal(auction_generate_data(mu = mu,
                                     alpha = alpha,
                                     sigma = sigma),
               NULL)
  # Requires numeric obs
  expect_equal(auction_generate_data(obs = "a",
                                     mu = mu,
                                     alpha = alpha,
                                     sigma = sigma),
              NULL)
  # Requires numeric n_bids
  expect_equal(auction_generate_data(obs = obs,
                                     mu = mu,
                                     alpha = alpha,
                                     sigma = sigma,
                                     n_bids = "a"),
               NULL)
  # N_bids sampled from 2:10
  data <- auction_generate_data(obs = obs,
                        mu = mu,
                        alpha = alpha,
                        sigma = sigma)
  expect_lte(min(data$n_bids), 10)
  expect_gte(min(data$n_bids), 2)
}
)

# Generating costs and proportional bid
test_that("Cost and bid generation",{
  obs = 200
  mu = 5
  alpha = 3
  sigma = .6
  n_bids = sample(2:10, obs, replace=TRUE)
  set.seed(2539)
  test_cost = auctionmodel:::auction__generate_cost(obs = obs,
                                        mu = mu,
                                        v.n = n_bids,
                                        alpha = alpha)
  # Numeric costs
  expect_equal(is.numeric(test_cost), TRUE)
  # Non-negative costs
  expect_gte(min(test_cost),
               0)
  # Check reasonable cost distribution
  plot(density(test_cost))

  # Check some specific parameter values
  mu = 0
  test_cost = auctionmodel:::auction__generate_cost(obs = obs,
                                                    mu = mu,
                                                    v.n = n_bids,
                                                    alpha = alpha)
  expect_equal(max(test_cost),
               0)
  expect_equal(min(test_cost),
               max(test_cost))

  mu = 5
  alpha = 1e10
  test_cost = auctionmodel:::auction__generate_cost(obs = obs,
                                                    mu = mu,
                                                    v.n = n_bids,
                                                    alpha = alpha)
  expect_equal(max(test_cost),
               5)
  expect_equal(min(test_cost),
               max(test_cost))

  # Non-negative bids
  gamma_1p1oa = gamma(1 + 1/alpha)
  test_w_bid = auctionmodel:::auction__generate_w_bid(obs = obs,
                                                      v.w_cost = test_cost,
                                                      v.n = n_bids,
                                                      mu = mu,
                                                      alpha = alpha,
                                                      gamma_1p1oa = gamma_1p1oa)
  # Numeric bids
  expect_equal(is.numeric(test_w_bid), TRUE)
  # Bids mark-up
  for(i in 1:obs){
    expect_gte(test_w_bid[i] - test_cost[i], 0)
  }
}
)

test_that("UH distribution",{
  obs = 100
  mu = 5
  alpha = 3
  sigma = .6

  test_u = auctionmodel:::auction__generate_u(obs = obs,
                                     sigma = sigma,
                                     u_dist = "dlnorm")
  # Numeric u
  expect_equal(is.numeric(test_u),
               TRUE)

  # Positive u
  expect_gte(min(test_u),
             0)

  # Log normal
  # Check similar distribution plots
  plot(density(test_u))
  curve(dlnorm(x, meanlog = -1/2 * log(1 + sigma^2), sdlog = sqrt(log(1 + sigma^2))), add = TRUE)

  # Specific parameter values
  sigma = 0
  test_u = auctionmodel:::auction__generate_u(obs = obs,
                                              sigma = sigma,
                                              u_dist = "dlnorm")
  expect_equal(min(test_u),
               1)
  expect_equal(max(test_u),
               min(test_u))


  # Gamma
  obs = 500
  mu = 5
  alpha = 3
  sigma = .6

  test_u = auctionmodel:::auction__generate_u(obs = obs,
                                              sigma = sigma,
                                              u_dist = "dgamma")
  # Check similar distribution plots
  plot(density(test_u))
  curve(dgamma(x, shape = 1/sigma^2, rate = 1/sigma^2), add = TRUE)

  # Specific parameter values
  sigma = 1
  set.seed(2547)
  test_u = auctionmodel:::auction__generate_u(obs = obs,
                                              sigma = sigma,
                                              u_dist = "dlnorm")
  gamma_test = rexp(obs, rate = 1)
  gamma_dist_test <- ks.test(test_u, gamma_test)
  expect_lte(gamma_dist_test$p.value,
             0.05)

  #Weibull
  obs = 500
  mu = 5
  alpha = 3
  sigma = .6
  set.seed(4892)
  test_u = auctionmodel:::auction__generate_u(obs = obs,
                                              sigma = sigma,
                                              u_dist = "dweibull")

  shape_solver <- stats::optimize(
    f = function(shape, sigma){
      return( abs( gamma(1+2/shape) / gamma(1+1/shape)^2 - 1 - sigma^2 ) )
    },
    interval = c(0.001, 127),
    tol = 1e-5,
    sigma = sigma
  )
  shape = shape_solver$minimum
  scale = 1/gamma(1+1/shape_solver$minimum)

  # Check similar distribution plots
  plot(density(test_u))
  curve(dweibull(x, shape = shape, scale = scale), add = TRUE)

  # Specific parameter values
  sigma = 0
  set.seed(2547)
  test_u = auctionmodel:::auction__generate_u(obs = 1e5,
                                              sigma = sigma,
                                              u_dist = "dweibull")
  expect_gte(min(test_u),
               0.9)
  expect_lte(max(test_u),
             1.1)

})

test_that("Winning bid",{
  obs = 500
  mu = 5
  alpha = 3
  sigma = .6
  test_u = auctionmodel:::auction__generate_u(obs = obs,
                                              sigma = sigma)
  n_bids = sample(2:10, obs, replace=TRUE)
  test_cost = auctionmodel:::auction__generate_cost(obs = obs,
                                                    mu = mu,
                                                    alpha = alpha,
                                                    v.n = n_bids)

  test_w_bid = auctionmodel:::auction__generate_w_bid(obs = obs,
                                                      v.w_cost = test_cost,
                                                      v.n = n_bids,
                                                      mu = mu,
                                                      alpha = alpha,
                                                      gamma_1p1oa = gamma(1+1/alpha))
  winning_bid = auctionmodel:::auction__generate_winning(all_x_vars = as.data.frame(rep(1, obs)),
                                                         beta = 0,
                                                         v.w_bid = test_w_bid,
                                                         v.u = test_u)

  expect_equal(winning_bid, test_w_bid*test_u)
})



