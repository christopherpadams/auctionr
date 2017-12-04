#######
#
# Examples on code use

library(auctionmodeling)

#############################
# Test default distribution
auction_v3(dat=auction__generate_data(obs=20), winning_bid = 'price', number_of_bids = 'num', num_cores = 3)

# Test weibull
auction_v3(dat=auction__generate_data(obs=20), winning_bid = 'price', number_of_bids = 'num', num_cores = 3,
           common_distributions = 'dweibull')

# Test with two distributions
res = auction_v3(common_distributions = c('dlnorm', 'dgamma'),
                 dat=auction__generate_data(obs=20),
                 winning_bid = 'price', number_of_bids = 'num',
                 num_cores = 3)
print(res)

# Test with two distributions
res = auction_v3(common_distributions = c('dlnorm', 'dgamma', 'dweibull'),
                 dat=auction__generate_data(obs=100),
                 winning_bid = 'price', number_of_bids = 'num',
                 num_cores = 3)
print(res)
