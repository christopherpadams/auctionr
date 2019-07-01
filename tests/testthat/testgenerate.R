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


