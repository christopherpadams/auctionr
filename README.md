Auction Modeling
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# auctionmodel

## Overview

A package for R to estimate private-value auction models while allowing
for unobservable auction-specific heterogeneity.

## Installation

``` r
# Install the development version from GitHub:

# install.packages("remotes")
remotes::install_github("alexjmac/auction_model")
```

## Getting started

There are two functions available in the package:

  - `auction_generate_data()` allows the user to generate sample data
    from the principal model used in the package.

  - `auction_model()` calculates maximum likelihood estimates of
    parameters of the principal model for the data provided by the user.

<!-- end list -->

``` r
library(auctionmodel)
set.seed(100)

dat <- auction_generate_data(obs = 100, mu = 10, alpha = 2, sigma = 0.2,
                             beta = c(-1,1), new_x_mean= c(-1,1), new_x_sd = c(0.5,0.8))

auction_model(dat,
              init_param =  c(8, 2, .5, .4, .6),
              num_cores = 1,
              method = "BFGS",
              control = list(trace=1, parscale = c(1,0.1,0.1,1,1)))
```

    ## initial  value 1339.327262 
    ## iter  10 value 434.301377
    ## iter  20 value 410.711195
    ## final  value 410.710822 
    ## converged

    ## $par
    ## [1] 11.0126727  1.7527689  0.2042297 -0.9206175  1.0680962
    ## 
    ## $value
    ## [1] 410.7108
    ## 
    ## $counts
    ## function gradient 
    ##       51       22 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL

## For further information

## License
