obs = 100000
alpha = 3
mu = 2
sigma = .1*mu + .1*alpha


data <- auction_generate_data(obs = obs,
                              mu = mu,
                              alpha = alpha,
                              n_bid = rep(5, obs),
                              sigma = sigma,
                              u_dist = "dlnorm")

data2 <- auction_generate_data(obs = obs,
                               mu = mu,
                               alpha = alpha,
                               n_bid = rep(5, obs),
                               sigma = sigma,
                               u_dist = "dweibull")
data3 <- auction_generate_data(obs = obs,
                               mu = mu,
                               alpha = alpha,
                               n_bid = rep(5, obs),
                               sigma = sigma,
                               u_dist = "dgamma")

plot(density(data$winning_bid), xlim = c(0, 5), ylim = c(0,1))
lines(density(data2$winning_bid), col = "red")
lines(density(data3$winning_bid), col = "blue")


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
curve(dweibull(x, shape = shape, scale = scale),
      xlim = c(0, 5), ylim = c(0,1), col = "red")
curve(dgamma(x, shape = 1/sigma^2, rate = 1/sigma^2),
      add = TRUE, col = "blue")
curve(dlnorm(x, meanlog = -1/2 * log(1 + sigma^2), sdlog = sqrt(log(1 + sigma^2))),
      add = TRUE)


