#######
#
# Examples on code use

library(auctionmodeling)

#############################
# Test default distribution
data <- auction_generate_data(obs=20)
model_auction(dat = data, winning_bid = 'price', number_of_bids = 'num', num_cores = 3)

# Test weibull
model_auction(dat=auction_generate_data(obs=20), winning_bid = 'price', number_of_bids = 'num', num_cores = 3,
           common_distributions = 'dweibull')

# Test with two distributions
res = model_auction(common_distributions = c('dlnorm', 'dgamma'),
                 dat=auction_generate_data(obs=20),
                 winning_bid = 'price', number_of_bids = 'num',
                 num_cores = 3)
print(res)

# Test with two distributions
res = model_auction(common_distributions = c('dlnorm', 'dgamma', 'dweibull'),
                 dat=auction_generate_data(obs=100),
                 winning_bid = 'price', number_of_bids = 'num',
                 num_cores = 3)
print(res)
