Auction Modeling
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# auctionr

## Overview

A package for R to estimate private-value auction models while allowing
for unobservable auction-specific heterogeneity.

## Installation

``` r
# Install a necessary dependency
install.packages("numDeriv")

# Download a ZIP archive of the repo from https://code.harvard.edu/HBS/rcs_amackay__auctionr
# Install the package
install.packages("[local path to the archive]", repos = NULL, type = "source")
```

## Getting started

There are two functions available in the package:

  - `auction_generate_data()` allows the user to generate sample data
    from the principal model used in the package.

  - `auction_model()` calculates maximum likelihood estimates of
    parameters of the principal model for the data provided by the user.

<!-- end list -->

``` r
library(auctionr)

set.seed(100)
dat <- auction_generate_data(obs = 100, mu = 10, alpha = 2, sigma = 0.2,
                             beta = c(-1,1), new_x_mean= c(-1,1), new_x_sd = c(0.5,0.8))

res <- auction_model(dat,
                    init_param =  c(8, 2, .5, .4, .6),
                    num_cores = 1,
                    method = "BFGS",
                    control = list(trace=1, parscale = c(1,0.1,0.1,1,1)),
                    std_err = TRUE)
```

    ## Running the optimizer using the BFGS method with starting values ( 8, 2, 0.5, 0.4, 0.6 )...
    ## 
    ## initial  value 1339.327262 
    ## iter  10 value 434.301377
    ## iter  20 value 410.711195
    ## final  value 410.710822 
    ## converged
    ## 
    ## Estimated parameters (SE):                              
    ##   mu      11.0127   (1.15264)  
    ##   alpha    1.75277  (0.185499) 
    ##   sigma    0.20423  (0.0352861)
    ##   beta[1] -0.920617 (0.0570396)
    ##   beta[2]  1.0681   (0.0400262)
    ## 
    ## Maximum log-likelihood = -410.711

## For further information

## License
