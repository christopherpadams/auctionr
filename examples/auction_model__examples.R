#' Example usage of the package 'auctionmodeling'
#'
#'   The code within will run the three main functions from the package:
#'
#'       auction_generate_data()
#'           generate sample data
#'
#'       auction_model()
#'           run full simulation on the data to find best-fitting distribution parameters
#'
#'       auction_model_likelihood()
#'           run only one iteration of the simulation, with specific distribution parameters
#'
#'
#'
#'
#' Comments
#' --------
#' Running a full simulation via function auction_model() may take considerable time.
#'
#'     The function auction_model() is used in example(s):
#'         Example 3: Simulation via function auction_model()
#'
#' It is recommended to first use function auction_model_likelihood() to verify
#' your data is correctly structured and that appropriate input arguments have
#' been set. This function will run only one iteration of the simulation, as
#' opposed to auction_model() which runs to convergence.
#'
#'     The function auction_model() is used in example(s):
#'         Example 2: Quick run via function auction_model_likelihood()

rm(list = ls())

library(auctionmodeling)



examples <- function(num_cores=4) {
  print('Running examples()')

  ### ------------------------------------------------------------
  ### Example 1: Data generation via function auction_generate_data()
  print('Example 1: Data generation via function auction_generate_data()')

  # Prepare source conditions
  #   Number of observations
  obs = 500

  #   Conditions for the Private Value distribution as well as Observed and Unobserved Heterogeneity
  mu = 5
  alpha = 2

  sigma = .5

  beta = c(.3, .2, .1, .4, .5)

  new_x_meanlog = c(2, 1, .5)
  new_x_sdlog = c(1, 1, 1)

  w = rlnorm(obs)
  x1 = rlnorm(obs) + .5*w
  x2 = .1*rlnorm(obs) + .3*w
  x_vars = data.frame(x1, x2)

  #   Number of bids
  v.n = sample(2:10, obs, replace=TRUE)

  # Generate source data
  print('  Generating data')
  data = auction_generate_data(obs = obs,
                               n_bids = v.n,
                               mu = mu,
                               alpha = alpha,
                               sigma = sigma,
                               beta = beta,
                               x_vars = x_vars,
                               new_x_meanlog = new_x_meanlog,
                               new_x_sdlog = new_x_sdlog)

  print('  summary(data)')
  print(summary(data))



  ### ------------------------------------------------------------
  ### Example 2: Quick run via function auction_model_likelihood()
  print('Example 2: Quick run via function auction_model_likelihood()')

  # Set input arguments related to the data
  #   Data
  #       Define the dataset and appropriate columns
  #         Use data from 'Example 1: Data generation...'
  #
  #   Private Value distribution
  #       Define parameters of the Private Value distribution
  #         Use same parameters as from 'Example 1: Data generation...'
  #
  #   Observed Heterogeneity
  #       Define 'beta' parameters of the Observed Heterogeneity
  #         Use same parameters as from 'Example 1: Data generation...'
  #
  #   Unobserved Heterogeneity
  #       Define distributions for unobserved heterogenity
  u_dist = c("dlnorm","dgamma")
  #       Define 'sigma' parameter of the Unobserved Heterogeneity
  #         Use same parameters as from 'Example 1: Data generation...'

  # Run function
  print('  Running auction_model_likelihood()')
  res = auction_model_likelihood(dat = data,
                                 winning_bid = "winning_bid",
                                 n_bids = "n_bids",
                                 u_dist = u_dist,
                                 alpha = alpha,
                                 sigma = sigma,
                                 beta = beta,
                                 num_cores = num_cores)

  print('  Result from "Running auction_model_likelihood()"')
  print(res)



  ### ------------------------------------------------------------
  ### Example 3: Simulation via function auction_model()
  print('Example 3: Simulation via function auction_model()')

  # Set input arguments related to the data
  #   Data
  #       Define the dataset and appropriate columns
  #         Use data from 'Example 1: Data generation...'
  #
  #   Private Value distribution
  #       Define the initial guesses for the parameters of the Private Value distribution
  #         We will skip this step and use package default initial guesses
  #
  #   Observed Heterogeneity
  #       Define the initial guesses for 'beta' parameters of the Observed Heterogeneity
  #         We will skip this step and use package default initial guesses
  #
  #   Unobserved Heterogeneity
  #       Define distributions for unobserved heterogenity
  u_dist = c("dlnorm","dgamma")
  #       Define the initial guess for the 'sigma' parameter of the Unobserved Heterogeneity
  #         We will skip this step and use package default initial guess

  # Set the non-data-related run arguments
  #   Report the status of the simulation every X iterations
  report = 100

  # Run function
  print('  Running auction_model()')
  res = auction_model(dat = data,
                      winning_bid = "winning_bid",
                      n_bids = "n_bids",
                      u_dist = u_dist,
                      num_cores = num_cores,
                      report = report)

  print('  Result from "Running auction_model()"')
  print(res)



  return(res)
}

# Run example code
res = examples(num_cores=4)
