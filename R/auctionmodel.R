#' Estimates a first-price auction model
#'
#'
#' @param dat data.frame containing the winning bids, number of bids, and \code{X} variables that describe the data.
#' @param winning_bid In list \code{dat}, the key whose value is a vector that holds the winning bids.
#' @param n_bids In list \code{dat}, the key whose value is a vector that holds the number of bids.
#' @param init_mu Value for \code{mu} for initial guess of the private value distribution.
#' @param init_alpha Value for \code{alpha} for initial guess of the private value distribution.
#' @param init_sigma Value for \code{sigma} for the initial guess of the unobserved heterogeneity.
#' @param init_beta Value for \code{beta} for initial guess of the private value distribution.
#' @param init_params Vector of init_mu, init_alpha, init_sigma, and init_beta, if not supplied separately
#' @param u_dist Which distributions to represent the unobserved heterogeneity.
#' @param num_cores The number of cores for running the model in parallel.
#' @param report Show optimization progress every X iterations, where X is defined by \code{report}
#'
#'
#' @details This function estimates a first-price auction model with conditional independent private values.
#' The model allows for unobserved heterogeneity that is common to all bidders in addition to observable
#' heterogeneity. The winning bid (Y) takes the form
#'
#' Y = B * U * h(X)
#'
#' where B Is the proportional winning bid, U is the unobserved heterogeneity, and h(X) controls for
#' observed heterogeneity. The model is log-linear so that
#' log(Y) = log(B) + log(U) + log(h(X)) and log(h(X)) = beta1 * X1 + beta2 * X2 + … .
#'
#' The (conditionally) independent private costs are drawn from a Weibull distribution
#' with parameters mu and alpha. The CDF of this distribution is given by
#'
#' F(c) = 1 – exp(- (c * 1/mu * Gamma(1 + 1/alpha))^(alpha))
#'
#' The unobserved heterogeneity can take on several different distributions, which are selected
#'  by the argument u_dist. It is normalized to have a mean of 1,
#'  with a free parameter sigma representing the standard deviation.
#'
#'
#' Representing the unobserved heterogeneity is controled by \code{u_dist}. This is either a string or vector of strings for
#' which distrubtions to use: \code{dlnorm} (default, if not supplied), \code{weibull}, and \code{gamma}.
#'
#' Either \code{ini_params} or the set \code{init_mu}, \code{init_alpha}, \code{init_beta}, and \code{init_sigma} must be supplied.
#' If \code{init_params} is supplied, the others are ignored.
#'
#' If more than one distrubtion is specified, the optimization is attempted for each distribution. Accordingly, the private value distribution
#' will change for each iteration of the optimization and again per distribution.
#'
#' The results are not displayed until all iterations and distributions are tested.
#'
#' This funtion utilizes the \code{Rsnow} framework within the \code{Rparallel} package. If \code{numcores} is not specified, this will be run using only
#' one CPU/core. One can use \code{parallel::detectCores()} to determine how many are available on your system, but you are not advised
#' to use all at once, as this may make your system unresponsive. Please see \code{Rparallel} and \code{Rsnow} for more details.
#'
#' VERIFY! Note that for supplied data, any rows with missing data will be removed prior to the numeric solver runs.
#'
#'
#' @return For each of the distributions speicifed in \code{u_dist}, ...
#'
#' @examples
#'
#' \dontrun{
#' DOES IT MAKE SENSE TO INCLUDE EXAMPLE OUTPUT??
#' PUT LONG EXAMPLE HERE
#' }
#'
#' @seealso \code{\link{auction_generate_data}}
#'
#'
#' @export

# Optimises f__ll_parallel
auction_model <- function(dat = NULL,
                          winning_bid = NULL,
                          n_bids = NULL,
                          init_mu = NULL,
                          init_alpha = NULL,
                          init_sigma = NULL,  #init_control
                          init_beta = NULL,   #init_common_sd
                          init_params = NULL,
                          u_dist = NULL, #common_distributions
                          num_cores = 1,
                          report=0
) {

  # Initialize environment
  hEnv_tmp = new.env()
  num_cores = auction__init_env(num_cores=num_cores, hEnv_tmp=hEnv_tmp)

  # Validate distributions requested for unobserved heterogeneity
  u_dist = auction__check__common_distrib(u_dist = u_dist)

  # Validate input data
  dat = auction__check_input_data(dat = dat,
                                  colName__winning_bid = winning_bid,
                                  colName__n_bids = n_bids)

  # Log the X terms within the input data
  dat[ ! names(dat) %in% c(winning_bid, n_bids) ] = log(
    dat[ ! names(dat) %in% c(winning_bid, n_bids) ] )

  # Prepare initial guesses
  vecInitGuess = auction__check_init_guess(dat = dat,
                                           init_mu = init_mu,
                                           init_alpha = init_alpha,
                                           init_sigma = init_sigma,
                                           init_beta = init_beta,
                                           init_params = init_params
  )

  # Prepare control parameters for numerical solver
  conv_ctrl = auction__get_conv_ctrl(vecInitGuess = vecInitGuess)

  # Set up parallelization of numerical solver
  cl = parallel::makeCluster(num_cores)

  parallel::clusterExport(cl,
                          varlist=c("vf__bid_function_fast",
                                    "vf__w_integrand_z_fast",
                                    "f__funk"),
                          envir = environment(auctionmodel) )

  # Run
  run_result = list()
  for (funcName in u_dist) {
    sFuncName = as.character(funcName)

    # Prepare tracker object
    hTracker = auction__tracker__build(hEnv_tmp=hEnv_tmp, report=report)

    # Run
    print(paste("Running |", sFuncName))
    #   Build function call parameter
    listFuncCall = list(funcName = sFuncName,
                        funcID = auction__get_id_distrib(sFuncName = sFuncName))
    #   Run
    run_result[[sFuncName]] = stats::optim(par=vecInitGuess,
                                           fn=f__ll_parallel,
                                           control=conv_ctrl,
                                           dat__winning_bid=dat[[winning_bid]],
                                           dat__n_bids=dat[[n_bids]],
                                           dat_X=t( dat[ ! names(dat) %in% c(winning_bid, n_bids) ] ),
                                           listFuncCall=listFuncCall,
                                           hTracker=hTracker,
                                           cl=cl)
    print(paste("End of run |", sFuncName))
  }
  # Release resources
  parallel::stopCluster(cl)
  # Prepare output
  res = auction__output_org(run_result=run_result,
                            dat_X__fields=names(dat)[ ! names(dat) %in% c(winning_bid, n_bids) ],
                            dat__winning_bid=dat[[winning_bid]])
  return(res)
}

#' NEW_TITLE - Estimates a first-price auction model
#'
#'
#' @param dat data.frame containing the winning bids, number of bids, and \code{X} variables that describe the data.
#' @param winning_bid In list \code{dat}, the key whose value is a vector that holds the winning bids.
#' @param n_bids In list \code{dat}, the key whose value is a vector that holds the number of bids.
#' @param mu Value for \code{mu} of the private value distribution.
#' @param alpha Value for \code{alpha} of the private value distribution.
#' @param sigma Value for \code{sigma} of the unobserved heterogeneity.
#' @param beta Value for \code{beta} of the private value distribution.
#' @param params Vector of mu, alpha, sigma, and beta, if not supplied separately
#' @param u_dist Which distributions to represent the unobserved heterogeneity.
#' @param num_cores The number of cores for running the model in parallel.
#'
#'
#' @details DETAILS_SECTION
#' #'
#' @return For each of the distributions specified in \code{u_dist}, ...
#'
#'
#' @export

# outputs log likelihood function
auction_model_likelihood <- function(dat = NULL,
                                     winning_bid = NULL,
                                     n_bids = NULL,
                                     mu = NULL,
                                     alpha = NULL,
                                     sigma = NULL,
                                     beta = NULL,
                                     params = NULL,
                                     u_dist = NULL,
                                     num_cores = 1 ) {

  # Initialize environment
  hEnv_tmp = new.env()
  num_cores = auction__init_env(num_cores=num_cores, hEnv_tmp=hEnv_tmp)

  # Validate distributions requested for unobserved heterogeneity
  u_dist = auction__check__common_distrib(u_dist = u_dist)

  # Validate input data
  dat = auction__check_input_data(dat = dat,
                                  colName__winning_bid = winning_bid,
                                  colName__n_bids = n_bids)

  # Log the X terms within the input data
  dat[ ! names(dat) %in% c(winning_bid, n_bids) ] = log(
    dat[ ! names(dat) %in% c(winning_bid, n_bids) ] )

  # Prepare initial guesses
  vecInitGuess = auction__check_init_guess(dat = dat,
                                           init_mu = mu,
                                           init_alpha = alpha,
                                           init_sigma = sigma,
                                           init_beta = beta,
                                           init_params = params
  )

  # Set up parallelization of numerical solver
  cl = parallel::makeCluster(num_cores)

  parallel::clusterExport(cl,
                          varlist=c("vf__bid_function_fast",
                                    "vf__w_integrand_z_fast",
                                    "f__funk"),
                          envir = environment(auctionmodel) )

  # Run
  run_result = list()
  for (funcName in u_dist) {
    sFuncName = as.character(funcName)

    # Prepare tracker object
    #   Set report=0 to ensure no reporting
    hTracker = auction__tracker__build(hEnv_tmp=hEnv_tmp, report=0)

    # Run
    print(paste("Running |", sFuncName))
    #   Build function call parameter
    listFuncCall = list(funcName = sFuncName,
                        funcID = auction__get_id_distrib(sFuncName = sFuncName))
    #   Run
    run_result[[sFuncName]] = f__ll_parallel(x0=vecInitGuess,
                                             dat__winning_bid=dat[[winning_bid]],
                                             dat__n_bids=dat[[n_bids]],
                                             dat_X=t( dat[ ! names(dat) %in% c(winning_bid, n_bids) ] ),
                                             listFuncCall=listFuncCall,
                                             hTracker=hTracker,
                                             cl=cl)
    print(paste("End of run |", sFuncName))
  }
  # Release resources
  parallel::stopCluster(cl)


  # Output from each run was negative log_likelihood
  #   Adjust
  for (keyResult in names(run_result)) {
    run_result[[keyResult]] = -run_result[[keyResult]]
  }

  return(run_result)
}

# output setup
auction__output_org <- function(run_result, dat_X__fields, dat__winning_bid) {

  # Initialize dataframe
  nCol = length(run_result)

  #   Initialize
  df = data.frame(matrix(ncol = length(run_result), nrow=0))
  #   Set columns
  colnames(df) = names(run_result)
  # Get indices for the 'par' data
  idxList = auction__x0_indices()
  # Fill
  for (funcName in names(run_result)) {
    sFuncName = as.character(funcName)

    # Private Values
    #   mu            VALUE
    #   a             VALUE
    df['Private Values', sFuncName] = ''
    df['  mu', sFuncName] = sprintf('%.4f',
                                    run_result[[sFuncName]]$par[idxList$pv_weibull_mu] )
    df['  a', sFuncName] = sprintf('%.4f',
                                   run_result[[sFuncName]]$par[idxList$pv_weibull_a] )

    # Unobserved Heterogeneity
    #   standard deviation        VALUE
    # Implied Parameters for U
    #   *  See manual *
    #   param 1                  [VALUE]
    #   param 2                  [VALUE]
    df['Unobserved Heterogeneity', sFuncName] = ''
    df['  standard deviation', sFuncName] = sprintf('%.4f', run_result[[sFuncName]]$par[idxList$unobs_dist_param] )
    df['Implied Parameters for U', sFuncName] = ''
    df['  * See manual *', sFuncName] = ''
    listParam = auction__get_unobs_params(distrib_std_dev =
                                            run_result[[sFuncName]]$par[idxList$unobs_dist_param],
                                          id_distrib =
                                            auction__get_id_distrib(
                                              sFuncName=sFuncName ) )
    for (iParam in 1:length(listParam)){
      # # df[sprintf('  param %d', iParam), sFuncName] = sprintf('%s = %.4f',
      # #                                                        names(listParam[iParam]),
      # #                                                        listParam[iParam] )
      # df[sprintf('  param %d', iParam), sFuncName] = sprintf('[%.4f (%s)]',
      #                                                        listParam[iParam],
      #                                                        names(listParam[iParam]) )
      df[sprintf('  param %d', iParam), sFuncName] = sprintf('[%.4f]',
                                                             listParam[iParam],
                                                             names(listParam[iParam]) )
    }

    # Observed Heterogeneity
    #   x1                    VALUE
    #   x2                    VALUE
    df['Observed Heterogeneity', sFuncName] = ''
    for (iX in 1:(1+length(run_result[[sFuncName]]$par)-idxList$x_terms__start)) {
      df[sprintf('  %s', dat_X__fields[iX]),
         sFuncName] = sprintf('%.4f',
                              run_result[[sFuncName]]$par[(idxList$x_terms__start+iX-1)] )
    }

    # Statistics
    #   log likelihood              VALUE
    #   Iterations                  VALUE
    #   StdDev: ln(Winning bids)    VALUE
    #   StdDev:                     VALUE
    df['Statistics', sFuncName] = ''
    df['  log likelihood', sFuncName] = sprintf('%.4f', run_result[[sFuncName]]$value )
    df['  Iterations', sFuncName] = run_result[[sFuncName]]$counts['function']
    df['  StdDev: ln(Winning bids)', sFuncName] = sprintf('%.4f', stats::sd(log(dat__winning_bid)))
    df['  StdDev: ln(Private Values)', sFuncName] = ''
    df['  StdDev: ln(Unobs. Hetero.)', sFuncName] = sprintf('%.4f', run_result[[sFuncName]]$par[idxList$unobs_dist_param])
    df['  StdDev: ln(Obs. Hetero.)', sFuncName] = ''
  }
  return(df)
}

#' Generate example data for running \code{\link{auctionmodel}}
#'
#'
#' @param obs Number of observations to draw
#'
#' @details This function generates example data for feeding into auctionmodel(). Specifically, the
#' winning bid, number of bids, and variables for the specified number of observations using random deviates of
#' the log normal distruction.
#'
#' @return A list with three keys:
#' \describe{
#' \item{price}{a numeric vector whose elements are the winning bids}
#' \item{num}{an integer vector containing the number of bids}
#' \item{x_terms}{a matrix of numeric elements of the x terms that describe the shape of the data}
#'}
#' For all three lists, the length of the key's values are the number of observations specified by the \code{obs} parameter.
#'
#' @examples
#' data <- auction_generate_data(100)
#' data$price
#' data$num
#' data$x_terms
#'
#' @seealso \code{\link{auctionmodel}}
#'
#'
#' @export
auction_generate_data <- function(obs = NULL,
                                  n_bids = NULL,
                                  mu = NULL,
                                  alpha = NULL,
                                  sigma = NULL,
                                  beta = NULL,
                                  params = NULL,
                                  u_dist = NULL,
                                  x_vars = NULL,
                                  new_x_meanlog = NULL,
                                  new_x_sdlog = NULL) {
  # Inspect parameters
  if (is.null(x_vars) && (is.null(new_x_meanlog))){
    beta = 0
  }
  if (is.null(params) && (is.null(mu)
                            || is.null(alpha)
                            || is.null(sigma)
                            || is.null(beta))) {
    print("Must specify either (mu, alpha, sigma, beta) or 'params'")
    return(NULL)
  } else if (is.null(obs)) {
    print("Must specify 'obs'")
    return(NULL)
  } else if ((! is.numeric(obs)) || (length(obs) != 1)) {
      print("'obs' must be numeric value")
      return(NULL)
  } else if (is.null(x_vars)){
    x_vars = as.data.frame(rep(1, obs))
    colnames(x_vars) = "x_vars"
  }
  # Inspect 'n_bids'
  #   'n_bids' may be a vector of number of bids.
  #   If specified, the length must be equal to
  #   obs. If not specified, they are drawn with
  #   replacement [sample(2:10, obs, replace=TRUE)]
  if (! is.null(n_bids)) {
    if ((! is.numeric(n_bids))
        || (! is.vector(n_bids))
        || (length(n_bids) != obs)) {
      print("'n_bids' must be vector, with length equal to 'obs")
      return(NULL)
    }
  } else {
    n_bids = sample(2:10, obs, replace=TRUE)
  }

  v.n = n_bids
  gamma_1p1oa = gamma(1 + 1/alpha)

  v.w_cost = auction__generate_cost(obs = obs,
                                    mu = mu,
                                    v.n = v.n,
                                    alpha = alpha)

  v.w_bid = auction__generate_w_bid(obs = obs,
                                    v.w_cost = v.w_cost,
                                    v.n = v.n,
                                    mu = mu,
                                    alpha = alpha,
                                    gamma_1p1oa = gamma_1p1oa)

  v.u = auction__generate_u(obs = obs, sigma = sigma)

  all_x_vars = auction__generate_x(obs = obs,
                                   x_vars = x_vars,
                                   beta = beta,
                                   new_x_meanlog = new_x_meanlog,
                                   new_x_sdlog = new_x_sdlog)

  v.winning_bid = auction__generate_winning(all_x_vars = all_x_vars,
                                            beta = beta,
                                            v.w_bid = v.w_bid,
                                            v.u = v.u)
  dat = data.frame(winning_bid = v.winning_bid, n_bids = v.n, all_x_vars)
  return(dat)

}

auction__generate_cost <- function(obs,
                                   mu,
                                   v.n,
                                   alpha) {
  # Winning cost
  v.w_cost = rep(NA, obs)
  # Generate cost
  for(i in 1:obs){
    costs = (mu/gamma(1+1/alpha))*(-log(1-stats::runif(v.n[i])))^(1/alpha)
    min_cost = min(costs)
    v.w_cost[i] = min(costs)
  }
  return(v.w_cost)
}

# v.w_bid[i] = f.bid_function(cost=min_cost, num_bids=v.n[i], mu=mu, alpha=alpha)
# Proportional bid function
auction__generate_w_bid <- function(obs,
                                    v.w_cost,
                                    v.n,
                                    mu,
                                    alpha,
                                    gamma_1p1oa) {
  v.w_bid = rep(NA, obs)
  for(i in 1:obs){
    v.w_bid[i] = f__bid_function_fast(cost=v.w_cost[i],
                                      n_bids=v.n[i],
                                      mu=mu,
                                      alpha=alpha,
                                      gamma_1p1oa=gamma_1p1oa)
  }
  return(v.w_bid)
}

# Unobserved heterogeneity
auction__generate_u <- function(obs,
                                sigma,
                                u_dist = "dlnorm") {
  for (funcName in u_dist) {
    sFuncName = as.character(funcName)
  }
  u_dist = auction__check__common_distrib(u_dist = u_dist)
  id_distrib = auction__get_id_distrib(sFuncName = sFuncName)
  listParam = auction__get_unobs_params(distrib_std_dev = sigma,
                                        id_distrib = id_distrib)
  if (id_distrib == 1) {
    #dgamma
    v.u = stats::rgamma(n = obs, shape = listParam$shape, rate = listParam$rate)
  }
  if (id_distrib == 2 || is.null(u_dist)) {
    v.u = stats::rlnorm(n = obs, meanlog = listParam$meanlog, sdlog = listParam$sdlog)
  }
  if (id_distrib == 3) {
    v.u = stats::rweibull(n = obs, shape = listParam$shape, scale = listParam$scale)
  }
  sdlog = sqrt(log(sigma^2 + 1))
  v.u = stats::rlnorm(n = obs, meanlog = -1/2*sdlog^2, sdlog = sdlog)
  return(v.u)
}

# Observed heterogeneity
auction__generate_x <- function(obs,
                                x_vars,
                                new_x_meanlog,
                                new_x_sdlog,
                                beta = beta) {
  # Inspect beta
  #   length of beta must be equal to
  #   the number of columns in x_vars
  #   plus the length of new_x_meanlog
    if ((! is.numeric(beta))
        || (! is.vector(beta))
        || (length(beta) !=
            ncol(x_vars) + length(new_x_meanlog))) {
      print(
         paste(
          "'beta' must be vector of length",
          "[# columns of 'x_vars'",
          "+ length of 'new_x_meanlog']"
        )
        )
      return(NULL)
    }

  # Inspect x_vars
  #   x_vars is a data.frame of control variables.
  #   The number of observations must be equal
  #   to obs.
  if ((!is.data.frame(x_vars))
      || (nrow(x_vars) != obs)) {
    print("'x_vars' must be dataframe with nrows='obs'")
    return(NULL)
  }

  # Inspect new_x_meanlog and new_x_sdlog
  if (!is.null(new_x_meanlog)) {
    if ((!is.numeric(new_x_meanlog))
        || (!is.vector(new_x_meanlog))) {
      print("'new_x_meanlog' must be numeric vector")
      return(NULL)
    } else if (is.null(new_x_sdlog)) {
      new_x_sdlog = rep(1, length(new_x_meanlog))
    } else if ((!is.numeric(new_x_sdlog))
               || (!is.vector(new_x_sdlog))
               ||
               (length(new_x_sdlog) != length(new_x_meanlog))) {
      print(paste(
        "'new_x_sdlog' must be numeric vector",
        "of same length as 'new_x_meanlog'"
      ))
      return(NULL)
    }
    # Generate new_x_vars
    new_x_vars = matrix(NA, obs, length(new_x_meanlog))
    colList = c()
    for (i.new_x in 1:length(new_x_meanlog)) {
      new_x_vars[, i.new_x] = stats::rlnorm(obs,
                                            meanlog = new_x_meanlog[i.new_x],
                                            sdlog = new_x_sdlog[i.new_x])
      colList = c(colList, paste0('X', i.new_x, '__rlnorm'))
    }
    colnames(new_x_vars) = colList
  } else if (is.null(new_x_meanlog)) {
    if (! is.null(new_x_sdlog)) {
      print("'new_x_sdlog' given but not 'new_x_meanlog'")
      return(NULL)
    }
    else {
      new_x_vars = rep(1, obs)
    }
  }

  # Gather all X terms
  all_x_vars = data.frame(x_vars, new_x_vars)
  return(all_x_vars)
}

# Calculate winning bid
auction__generate_winning <- function(all_x_vars,
                                      beta,
                                      v.w_bid,
                                      v.u) {
  v.h_x = exp(colSums(beta*t(log(all_x_vars))))
  v.winning_bid = v.w_bid*v.u*v.h_x
  return(v.winning_bid)
}



auction__gen_err_msg <- function(res) {
  # Goal: Print out an error message and then stops execution of the main script

  errMsg = paste('\n\tError Code=', res['err_code'], '\n\tError msg=', res['err_msg'], sep='')
  stop(errMsg)
}

auction__init_env <- function(num_cores, hEnv_tmp) {
  # Goal:
  #   - Load all required packages, stop execution is any packages are missing
  #   - Check number of cores requested

  # Load required packages
  auction__load_packages()
  # Check number of cores requested
  num_cores = auction__check__num_cores(num_cores = num_cores)

  # Create reference for tracking numeric-solver progress
  auction__tracker__create_class(hEnv_tmp)

  return(num_cores)
}

auction__load_packages <- function () {
  # Goal: Load all required packages, stop execution if any are missing

  # Which packages are required?
  listRequiredPackages = c('parallel')
  # Try to load each required package, record which ones fail to load
  listMissingPackages = c()
  for (reqPkg in listRequiredPackages) {
    if ( ! require(reqPkg, character.only=TRUE)) {
      listMissingPackages = c(listMissingPackages, reqPkg)
    }
  }
  # Report which packages failed to load, if any
  if ( length(listMissingPackages) > 0 ) {
    res = list()
    res['err_code'] = 1
    res['err_msg'] = paste0("Unable to load the following packages: ",
                            paste(listMissingPackages, collapse=','))

    auction__gen_err_msg(res)
  }
}

auction__check__num_cores <- function(num_cores) {
  # Goal: Check number of cores requested

  # If no # workers requested, then set to 1
  if (is.null(num_cores)) {
    print(paste0("Setting number of cores to 1"))
    num_cores = 1
  }
  # Check if valid input received
  if (
    ( ! is.numeric(num_cores) ) ||
    ( num_cores%%1 != 0) ||
    ( num_cores < 0 )
  ) {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'num_cores' | Must be natural number"
    auction__gen_err_msg(res)
  } else {
    # Check cores available
    num_cores__avail = parallel::detectCores()
    if ( num_cores > num_cores__avail ) {
      print(paste0("Warning: You have requested ", num_cores,
                   " cores but only have ", num_cores__avail,
                   " cores available"))

      # Adjust number of workers we will request
      num_cores = num_cores__avail
      print(paste0("\tSetting # parallel workers to", num_cores))
    }
  }
  return(num_cores)
}

auction__check__common_distrib <- function(u_dist) {
  if (is.null(u_dist)) {
    u_dist = 'dlnorm'
  } else if (is.character(u_dist)) {
    for (distrib in u_dist) {
      if ((distrib != 'dlnorm') && (distrib != 'dweibull') && (distrib != 'dgamma')) {
        res = list()
        res['err_code'] = 2
        res['err_msg'] = paste("Invalid input for 'u_dist' | Entry '", distrib, "' is invalid", sep='')
        auction__gen_err_msg(res)
      }
    }
  } else {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'u_dist' | Must be string or vector of strings"
    auction__gen_err_msg(res)
  }
  return(u_dist)
}

auction__check_input_data <- function(dat, colName__winning_bid, colName__n_bids) {
  # Goal: Make sure the input data has required columns

  # Ensure 'dat' is a dataframe
  if (! is.data.frame(dat)) {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = paste0("Unable to find the following columns within input data: ",
                            paste(listMissingColName, collapse=','))
    auction__gen_err_msg(res)
  } else {

    # Ensure 'dat' has all required columns
    #   Get list of required column names
    colList_req = c(colName__winning_bid, colName__n_bids)

    #   Get list of column names
    colList = colnames(dat)
    #   Compare list against list of required columns
    listMissingColName = c()
    for (reqCol in colList_req) {
      if ( ! reqCol %in% colList ) {
        listMissingColName = c(listMissingColName, reqCol)
      }
    }
    if ( length(listMissingColName) > 0 ) {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = paste0("Unable to find the following columns within input data: ",
                              paste(listMissingColName, collapse=','))
      auction__gen_err_msg(res)
    } else {

      # Remove non-numeric columns
      for (colName in colList) {
        if ((class(dat[[colName]])=='factor') || (mode(dat[[colName]]) != 'numeric')) {
          # Drop column
          #   If price or number of bids, then abort
          if ((colName == colName_price) || (colName == colName_num)){
            # Abort
            res = list()
            res['err_code'] = 2
            res['err_msg'] = "'winning bid' and/or 'number of bids' has non-numeric data"
            auction__gen_err_msg(res)
          } else {
            # Remove
            dat[colName] = NULL
            colList = colList[colName != colList]
          }

        }
      }

      # Remove rows with any missing data
      #   Identify rows
      listRow = c()
      if (is.data.frame(dat)) {
        # Dealing with a dataframe
        for (iRow in 1:dim(dat)[1]) {
          if (any(is.na(dat[iRow,]))) {
            listRow = c(listRow, iRow)
          }
        }
      }
      if (length(listRow) > 0) {
        # Remove these rows
        dat = dat[-listRow,]
      }

      return(dat)
    }
  }
  return(NULL)
}

auction__check_init_guess <- function(dat = dat,
                                      init_mu,
                                      init_alpha,
                                      init_sigma,
                                      init_beta,
                                      init_params) {

  # Find number of X parameters within dat
  #   (Total number of columns) - (column winning bid) - (column number of bids)
  nCol_X = auction__get_num_columns__dat(dat) - 2
  # Initialize initial guess vector
  x0 = numeric(2 + 1 + nCol_X)
  # Fill vector
  #   Gather default values
  def_pv_mu = 8
  def_pv_a = 2
  def_unobs_stddev = 0.5
  def_x = rep(0.5, nCol_X)
  #   Gather position indices
  idxList = auction__x0_indices()


  # Check to see if 'init_params' is given
  if (! is.null(init_params)) {
    print("init_params given, other initial guesses will be ignored")

    # Must be vector with length 3 + number of X variables
    if (! is.vector(init_params)) {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_params' | must be a vector"
      auction__gen_err_msg(res)
    } else if (length(init_params) != length(x0)) {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_params' | invalid vector length"
      auction__gen_err_msg(res)
    } else {
      x0 = init_params
    }
  } else {

    if (is.null(init_mu)) {
      x0[idxList$pv_weibull_mu] = def_pv_mu
    } else if (is.numeric(init_mu)) {
      x0[idxList$pv_weibull_mu] = init_mu
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_mu'"
      auction__gen_err_msg(res)
    }

    if (is.null(init_alpha)) {
      x0[idxList$pv_weibull_a] = def_pv_a
    } else if (is.numeric(init_alpha)) {
      x0[idxList$pv_weibull_a] = init_alpha
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_alpha'"
      auction__gen_err_msg(res)
    }

    # if (is.null(init_beta)) {
    #   x0[idxList$unobs_dist_param] = def_pv_mu
    # } else if (is.numeric(init_beta)) {
    #
    #   print(idxList$unobs_dist_param)
    #   print(init_beta)
    #
    #   x0[idxList$unobs_dist_param] = init_beta
    # } else {
    #   res = list()
    #   res['err_code'] = 2
    #   res['err_msg'] = "Invalid input for 'init_beta'"
    #   auction__gen_err_msg(res)
    # }
    #
    # if (is.null(init_sigma)) {
    #   x0[idxList$x_terms__start:length(x0)] = def_x
    # } else if (is.numeric(init_sigma)) {
    #   x0[idxList$x_terms__start:length(x0)] = init_sigma
    # } else {
    #   res = list()
    #   res['err_code'] = 2
    #   res['err_msg'] = "Invalid input for 'init_sigma'"
    #   auction__gen_err_msg(res)
    # }

    if (is.null(init_sigma)) {
      x0[idxList$unobs_dist_param] = def_pv_mu
    } else if (is.numeric(init_sigma)) {
      x0[idxList$unobs_dist_param] = init_sigma
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_sigma'"
      auction__gen_err_msg(res)
    }

    if (is.null(init_beta)) {
      x0[idxList$x_terms__start:length(x0)] = def_x
    } else if (is.numeric(init_beta)) {
      x0[idxList$x_terms__start:length(x0)] = init_beta
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_beta'"
      auction__gen_err_msg(res)
    }

  }
  return(x0)
}

auction__get_conv_ctrl <- function(vecInitGuess) {
  # Max number of iterations = maxit
  maxit = 2000

  # Step sizes
  #   Initialize
  parscale = numeric(length(vecInitGuess))
  #   Gather position indices
  idxList = auction__x0_indices()
  #   Gather default values
  def_pv_mu = 1
  def_pv_a = 0.1
  def_unobs_stddev = 1
  def_x = c(0.1, rep(1, (length(vecInitGuess)-idxList$x_terms__start)) )
  #
  #   Fill
  parscale[idxList$pv_weibull_mu] = def_pv_mu
  parscale[idxList$pv_weibull_a] = def_pv_a
  parscale[idxList$unobs_dist_param] = def_unobs_stddev
  parscale[idxList$x_terms__start:length(vecInitGuess)] = def_x

  return( list(maxit = maxit, parscale = parscale ) )
}

auction__get_id_distrib <- function(sFuncName) {
  if (sFuncName == 'dgamma') {
    id_distrib = 1
  } else if (sFuncName == 'dlnorm') {
    id_distrib = 2
  } else if (sFuncName == 'dweibull') {
    id_distrib = 3
  } else {
    res = list()
    res['err_code'] = 3
    res['err_msg'] = paste("Unable to get id_distrib for '", sFuncName, "'", sep='')
    auction__gen_err_msg(res)
  }
  return(id_distrib)
}

auction__get_unobs_params <- function(distrib_std_dev, id_distrib) {
  if (id_distrib == 1) {
    # dgamma
    listParam = auction__get_distrib_params__gamma(distrib_std_dev = distrib_std_dev)
  } else if (id_distrib == 2) {
    # dlnorm
    listParam = auction__get_distrib_params__lognorm(distrib_std_dev = distrib_std_dev)
  } else if (id_distrib == 3) {
    # dweibull
    listParam = auction__get_distrib_params__weibull(distrib_std_dev = distrib_std_dev)
  } else {
    res = list()
    res['err_code'] = 3
    res['err_msg'] = paste("Unknown id_distrib '", id_distrib, "'", sep='')
    auction__gen_err_msg(res)
  }
  return(listParam)
}

auction__get_distrib_params__lognorm <- function(distrib_std_dev) {
  # Given std dev and E(X) = 1, calculate meanlog and sdlog
  tmp = log(1+distrib_std_dev^2)
  return(list(sdlog=sqrt(tmp), meanlog=-1/2*tmp))
  # return(list( meanlog=(-distrib_std_dev^2*1/2), sdlog = distrib_std_dev ))
}

auction__get_distrib_params__weibull <- function(distrib_std_dev) {
  # Given std dev and E(X) = 1, calculate scale and shape
  #   S^2 + 1 = GAMMA(1+2/shape) / [GAMMA(1+1/shape)]^2
  #     need to numerically solve
  res_solver = stats::optimize(
    f = function(shape, S){
      return( abs( gamma(1+2/shape) / gamma(1+1/shape)^2 - 1 - S^2 ) )
    },
    interval = c(0.001, 127),
    tol = 1e-5,
    S = distrib_std_dev
  )
  return(list(shape = res_solver$minimum, scale = 1/gamma(1+1/res_solver$minimum)))
}

auction__get_distrib_params__gamma <- function(distrib_std_dev) {
  # Given std dev and E(X) = 1, calculate rate and shape
  tmp = 1/distrib_std_dev^2
  return(list(shape = tmp, rate = tmp))
}

auction__get_private_value_stats <- function(weibull_scale, weibull_shape) {
  std_dev = sqrt(
    weibull_scale^2 * (
      gamma(1 + 2 / weibull_shape) - (gamma(1 + 1 / weibull_shape))^2
    )
  )

  return(std_dev)
}

auction__get_num_columns__dat <- function(dat) {
  if (is.data.frame(dat)) {
    nParams_dat = dim(dat)[2]
  } else {
    # Expected to be list, with each key-value pairing potentially having more than 1 column
    #   e.g. either "x1", "x2", ... or "x_terms"
    nParams_dat = 0
    for (colName in names(dat)) {
      nCol = dim(dat[[colName]])[2]
      if (is.null(nCol)) {
        nCol = 1
      }
      nParams_dat = nParams_dat + nCol
    }
  }
  return(nParams_dat)
}

auction__x0_indices <- function() {
  return( list(
    pv_weibull_mu = 1,
    pv_weibull_a = 2,
    unobs_dist_param = 3,
    x_terms__start = 4
  ) )
}

auction__tracker__create_class <-function(hEnv_tmp) {

  # with(hEnv_tmp, {
  #   setRefClass("auctionmodel__tracker",
  #               fields=list(
  #                 iter="numeric",
  #                 report="numeric"
  #               ))
  # } )

  hEnv_tmp$iter = 0
  hEnv_tmp$report = 0

}

auction__tracker__build <-function (hEnv_tmp, report) {
  # # 'report' must be a numeric, greater than or equal to 0
  # #   round 'report' just in case
  # if ((is.null(report)) || (! is.numeric(report)) || report < 0) {
  #   report = 0
  # } else {
  #   report = round(report)
  # }
  # # build
  # hTracker = new("auctionmodel__tracker", iter=0, report=report)

  hEnv_tmp$iter = 0
  hEnv_tmp$report = report

  # return(hTracker)
  return(hEnv_tmp)
}

# log likelihood of winning bids
f__ll_parallel = function(x0, dat__winning_bid, dat__n_bids, dat_X, listFuncCall, hTracker, cl){

  listIdx = auction__x0_indices()

  if(x0[listIdx$unobs_dist_param] <= 0.1) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_mu] <= 0) > 0) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_a] <= 0.01) > 0) return(-Inf)

  param.h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1)]
  v.h = exp( colSums(param.h * dat_X) )

  listFuncCall$argList = auction__get_unobs_params(
    distrib_std_dev = x0[listIdx$unobs_dist_param],
    id_distrib = listFuncCall$funcID)

  # calculate bid density
  v.f_w = parallel::parApply(cl = cl,
                             X = cbind(dat__winning_bid/v.h,
                                       dat__n_bids,
                                       x0[listIdx$pv_weibull_mu],
                                       x0[listIdx$pv_weibull_a],
                                       gamma(1 + 1/x0[listIdx$pv_weibull_a]) ),
                             MARGIN = 1,
                             FUN = f__funk,
                             listFuncCall=listFuncCall)
  # density of winning bid (w/ controls)
  v.f_y = v.f_w/v.h
  # log likelihood for observed (winning) bid
  log_likelihood = -sum(log(v.f_y))

  if (hTracker$report != 0) {
    hTracker$iter = hTracker$iter + 1
    if (hTracker$iter %% hTracker$report == 0) {
      # Print to user
      #   Iteration
      #   log-likelihood [2-digits]
      #   parameter values that are passed to f__ll_parallel() [4-digits]
      #     x0
      #       mu, alpha, sd_U, beta(s)
      #     listFuncCall
      #       function name

      sDebug = sprintf('Iteration %d | distrib %s | log-likelihood %.2f | mu %.4f | alpha %.4f | sd_U %.4f | beta(s) %s',
                       hTracker$iter,
                       listFuncCall$funcName,
                       log_likelihood,
                       x0[listIdx$pv_weibull_mu],
                       x0[listIdx$pv_weibull_a],
                       x0[listIdx$unobs_dist_param],
                       paste(
                         sapply(listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1),
                                function(i) sprintf('%.4f', x0[i])),
                         sep='', collapse=', '))
      print(sDebug)
    }
  }

  return(log_likelihood)
}

# bid density (integrate over vf__w_integrand_z_fast)
f__funk = function(data_vec, listFuncCall) {
  val = stats::integrate(vf__w_integrand_z_fast, w_bid=data_vec[1],
                         n_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                         gamma_1p1oa=data_vec[5], listFuncCall=listFuncCall, lower=0, upper=Inf,
                         abs.tol = 1e-10)
  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

# integrand for density of observed bid, given proportional bid function, cost dist., UH dist.
vf__w_integrand_z_fast = function(z, w_bid, n_bids, mu, alpha, gamma_1p1oa, listFuncCall){

  b_z = vf__bid_function_fast(cost=z, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa)

  listFuncCall$argList$x = w_bid/b_z

  vals = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    do.call(
      match.fun(listFuncCall$funcName),
      listFuncCall$argList
    )

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-n_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}

# proportional bid function given cost, pv dist parameters, no bids
f__bid_function_fast = function(cost, n_bids, mu, alpha, gamma_1p1oa){

  if (exp(-(n_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha) == 0) {
    return(cost + mu/alpha*(n_bids-1)^(-1/alpha)*1/gamma_1p1oa*
             ((n_bids-1)*(gamma_1p1oa/mu*cost)^alpha)^(1/alpha-1))
  }

  cost + 1/alpha*(mu/gamma_1p1oa)*(n_bids-1)^(-1/alpha)*
    stats::pgamma((n_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(n_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha)
}

vf__bid_function_fast = Vectorize(FUN = f__bid_function_fast,vectorize.args = "cost")
