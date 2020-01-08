context("Auction model")

library(auctionmodel)

test_that("Requires appropriate parameters", {
  obs = 100
  mu = 5
  alpha = 1
  sigma = .6
  exd <- readRDS("exd.rds")

  expect_error(auction_model(exd),
               "Argument.+required")

  # Requires numeric inputs
  expect_error(auction_model(exd, init_param = c(mu, alpha, sigma, 0, 0, 0, 0, 0)),
               "column.+must be numeric")

  # Must have init_param values for everything
  expect_error(auction_model(exd, init_param = c(mu, alpha, sigma, 0, 0, 0, 0)),
               "Argument.+must be of length")

  # This one should work
  expect_error(m1 <- auction_model(exd[, 1:5], init_param = c(mu, alpha, sigma, 0, 0, 0)),
               NA)
  expect_equal(length(m1$par), 6)
}
)
