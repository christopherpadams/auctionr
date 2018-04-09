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
#' @param u_distributions Which distributions to represent the unobserved heterogeneity.
#' @param num_cores The number of cores for running the model in parallel.
#' @param debug
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
#' Representing the unobserved heterogeneity is controled by \code{u_distributions}. This is either a string or vector of strings for
#' which distrubtions to use: \code{dlnorm} (default, if not supplied), \code{weibull}, and \code{weibull}.
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
#' @return For each of the distributions speicifed in \code{u_distributions}, ...
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
auction_model <- function(dat = NULL,
                       winning_bid = NULL,
                       n_bids = NULL,
                       init_mu = NULL,
                       init_alpha = NULL,
                       init_sigma = NULL,  #init_control
                       init_beta = NULL,   #init_common_sd
                       init_params = NULL,
                       u_distributions = NULL, #common_distributions
                       num_cores = 1
                       ) {
  # Initialize environment
  num_cores = auction__init_env(num_cores=num_cores)

  # Validate distributions requested for unobserved heterogeneity
  u_distributions = auction__check__common_distrib(u_distributions = u_distributions)

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
                envir = environment(auction_model)
                )

  # Run
  run_result = list()
  for (funcName in u_distributions) {
    sFuncName = as.character(funcName)

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
                                           cl=cl)
    print(paste("End of run |", sFuncName))
  }
  # Release resources
  parallel::stopCluster(cl)
  # Prepare output
  res = auction__output_org(run_result=run_result,
                            dat_X__fields = names(dat)[ ! names(dat) %in% c(winning_bid, n_bids) ])
  return(res)
}

auction__output_org <- function(run_result, dat_X__fields) {

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

    df['Private_Value', sFuncName] = ''
    df['  mu', sFuncName] = sprintf('%.4f',
                                    run_result[[sFuncName]]$par[idxList$pv_weibull_mu] )
    df['  a', sFuncName] = sprintf('%.4f',
                                   run_result[[sFuncName]]$par[idxList$pv_weibull_a] )

    df['Unobserved_Hetero', sFuncName] = ''
    listParam = auction__get_unobs_params(distrib_std_dev =
                                            run_result[[sFuncName]]$par[idxList$unobs_dist_param],
                                          id_distrib =
                                            auction__get_id_distrib(
                                              sFuncName=sFuncName ) )
    for (iParam in 1:length(listParam)){
      df[sprintf('  param %d', iParam), sFuncName] = sprintf('%s = %.4f',
                                                             names(listParam[iParam]),
                                                             listParam[iParam] )
    }

    df['X_terms', sFuncName] = ''
    for (iX in 1:(1+length(run_result[[sFuncName]]$par)-idxList$x_terms__start)) {
      df[sprintf('  %s', dat_X__fields[iX]),
         sFuncName] = sprintf('%.5f',
                              run_result[[sFuncName]]$par[(idxList$x_terms__start+iX-1)] )
    }


    df['Statistics', sFuncName] = ''
    df['  inv__log_likelihood', sFuncName] = sprintf('%.5f',
                                                     run_result[[sFuncName]]$value )
    df['  StdDev__Private_Value', sFuncName] = sprintf('%.5f',
                                                       auction__get_private_value_stats(
                                                         weibull_scale=run_result[[sFuncName]]$par[idxList$pv_weibull_a],
                                                         weibull_shape=run_result[[sFuncName]]$par[idxList$pv_weibull_mu]
                                                       ))
    df['  StdDev__unobserved_hetero', sFuncName] = sprintf('%.5f',
                                                           run_result[[sFuncName]]$par[idxList$unobs_dist_param] )
    df['  Iterations', sFuncName] = run_result[[sFuncName]]$counts['function']

  }
  return(df)
}

#' Generate example data for running \code{\link{auction_model}}
#'
#'
#' @param obs Number of observations to draw
#'
#' @details This function generates example data for feeding into auction_model(). Specifically, the
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
#' @seealso \code{\link{auction_model}}
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
  if (is.null(params) && (is.null(mu)
                          || is.null(alpha)
                          || is.null(sigma)
                          || is.null(beta))) {
    print("Must specify either (mu, alpha, sigma, beta) or 'params'")
    return(NULL)
  } else if (is.null(obs)) {
    print("Must specify 'obs'")
    return(NULL)
  } else {

    if ((! is.numeric(obs)) || (length(obs) != 1)) {
      print("'obs' must be numeric value")
      return(NULL)
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

    # Inspect x_vars
    #   x_vars is a data.frame of control variables.
    #   The number of observations must be equal
    #   to obs.
    if ((! is.data.frame(x_vars))
        || (nrow(x_vars) != obs)) {
      print("'x_vars' must be dataframe with nrows='obs'")
      return(NULL)
    }

    # Inspect beta
    #   length of beta must be equal to
    #   the number of columns in x_vars
    #   plus the length of new_x_meanlog
    if ((! is.numeric(beta))
        || (! is.vector(beta))
        || (length(beta) !=
            ncol(x_vars) + length(new_x_meanlog))) {
      print(paste("'beta' must be vector of length",
                  "[# columns of 'x_vars'",
                  "+ length of 'new_x_meanlog']"))
      return(NULL)
    }

    # Inspect new_x_meanlog and new_x_sdlog
    if (! is.null(new_x_meanlog)) {
      if ((! is.numeric(new_x_meanlog))
          || (! is.vector(new_x_meanlog))) {
        print("'new_x_meanlog' must be numeric vector")
        return(NULL)
      } else if (is.null(new_x_sdlog)) {
        new_x_sdlog = rep(1, length(new_x_meanlog))
      } else if ((! is.numeric(new_x_sdlog))
                 || (! is.vector(new_x_sdlog))
                 || (length(new_x_sdlog) != length(new_x_meanlog))) {
        print(paste("'new_x_sdlog' must be numeric vector",
                    "of same length as 'new_x_meanlog'"))
        return(NULL)
      }
    } else {
      if (! is.null(new_x_sdlog)) {
        print("'new_x_sdlog' given but not 'new_x_meanlog'")
        return(NULL)
      }
    }

  }


  # Number of bids
  v.n = n_bids

  # Winning bids
  v.w_bid = rep(NA, obs)
  v.w_cost = rep(NA, obs)
  gamma_1p1oa = gamma(1 + 1/alpha)
  for(i in 1:obs){
    costs = (mu/gamma(1+1/alpha))*(-log(1-runif(v.n[i])))^(1/alpha)
    min_cost = min(costs)
    v.w_cost[i] = min(costs)
    # v.w_bid[i] = f.bid_function(cost=min_cost, num_bids=v.n[i], mu=mu, alpha=alpha)
    v.w_bid[i] = f__bid_function_fast(winning_bid=min_cost,
                                      n_bids=v.n[i],
                                      mu=mu,
                                      alpha=alpha,
                                      gamma_1p1oa=gamma_1p1oa)
  }
  stopifnot(mean(is.na(v.w_bid)) == 0)




  # Unobserved heterogeneity
  sdlog = sqrt(log(sigma^2+1))
  v.u = rlnorm(n = obs, meanlog = -1/2*sdlog^2, sdlog = sdlog )

  # Observed heterogeneity
  new_x_vars = matrix(NA, obs, length(new_x_meanlog))
  colList = c()
  for(i.new_x in 1:length(new_x_meanlog)){
    new_x_vars[, i.new_x] = rlnorm(obs,
                                   meanlog = new_x_meanlog[i.new_x],
                                   sdlog = new_x_sdlog[i.new_x])
    colList = c(colList, paste0('X', i.new_x, '__rlnorm'))
  }
  colnames(new_x_vars) = colList

  # Gather all X terms
  all_x_vars = data.frame(x_vars, new_x_vars)
  # Calculate winning bid
  v.h_x = exp(colSums(beta*t(log(all_x_vars))))
  v.winning_bid = v.w_bid*v.u*v.h_x
  # Build dataframe
  dat = data.frame(winning_bid = v.winning_bid, n_bids = v.n, all_x_vars)
  return(dat)
}

auction__gen_err_msg <- function(res) {
  # Goal: Print out an error message and then stops execution of the main script

  errMsg = paste('\n\tError Code=', res['err_code'], '\n\tError msg=', res['err_msg'], sep='')
  stop(errMsg)
}

auction__init_env <- function(num_cores) {
  # Goal:
  #   - Load all required packages, stop execution is any packages are missing
  #   - Check number of cores requested

  # Load required packages
  auction__load_packages()
  # Check number of cores requested
  num_cores = auction__check__num_cores(num_cores = num_cores)
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

auction__check__common_distrib <- function(u_distributions) {
  if (is.null(u_distributions)) {
    u_distributions = 'dlnorm'
  } else if (is.character(u_distributions)) {
    for (distrib in u_distributions) {
      if ((distrib != 'dlnorm') && (distrib != 'dweibull') && (distrib != 'dgamma')) {
        res = list()
        res['err_code'] = 2
        res['err_msg'] = paste("Invalid input for 'u_distributions' | Entry '", distrib, "' is invalid", sep='')
        auction__gen_err_msg(res)
      }
    }
  } else {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'u_distributions' | Must be string or vector of strings"
    auction__gen_err_msg(res)
  }
  return(u_distributions)
}

auction__check_input_data <- function(dat, colName__winning_bid, colName__n_bids) {
  # Goal: Make sure the input data has required columns

  # Ensure 'dat' is a dataframe
  if (! is.data.frame(dat)) {
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

    if (is.null(init_beta)) {
      x0[idxList$unobs_dist_param] = def_pv_mu
    } else if (is.numeric(init_beta)) {
      x0[idxList$unobs_dist_param] = init_beta
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_beta'"
      auction__gen_err_msg(res)
    }

    if (is.null(init_sigma)) {
      x0[idxList$x_terms__start:length(x0)] = def_x
    } else if (is.numeric(init_sigma)) {
      x0[idxList$x_terms__start:length(x0)] = init_sigma
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_sigma'"
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

f__ll_parallel = function(x0, dat__winning_bid, dat__n_bids, dat_X, listFuncCall, cl){
  listIdx = auction__x0_indices()

  if(x0[listIdx$unobs_dist_param] <= 0.1) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_mu] <= 0) > 0) return(-Inf)
  if(sum(x0[listIdx$pv_weibull_a] <= 0.01) > 0) return(-Inf)

  param.h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[1]-1)]
  v.h = exp(colSums(param.h*dat_X))

  listFuncCall$argList = auction__get_unobs_params(
    distrib_std_dev = x0[listIdx$unobs_dist_param],
    id_distrib = listFuncCall$funcID)

  v.f_w = parallel::parApply(cl = cl,
                             X = cbind(dat__winning_bid/v.h,
                                       dat__n_bids,
                                       x0[listIdx$pv_weibull_mu],
                                       x0[listIdx$pv_weibull_a],
                                       gamma(1 + 1/x0[listIdx$pv_weibull_a]) ),
                             MARGIN = 1,
                             FUN = f__funk,
                             listFuncCall=listFuncCall)
  v.f_y = v.f_w/v.h
  return(-sum(log(v.f_y)))
}

f__funk = function(data_vec, listFuncCall){
  val = integrate(vf__w_integrand_z_fast, w_bid=data_vec[1],
                  n_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], listFuncCall=listFuncCall, lower=0, upper=Inf,
                  abs.tol = 1e-10)
  if(val$message != "OK") stop("Integration failed.")
  return(val$value)
}

vf__w_integrand_z_fast = function(z, w_bid, n_bids, mu, alpha, gamma_1p1oa, listFuncCall){

  b_z = vf__bid_function_fast(winning_bid=z, n_bids=n_bids, mu=mu, alpha=alpha, gamma_1p1oa)

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

f__bid_function_fast = function(winning_bid, n_bids, mu, alpha, gamma_1p1oa){

  if (exp(-(n_bids-1)*(1/(mu/gamma_1p1oa)*winning_bid)^alpha) == 0) {
    return(winning_bid + mu/alpha*(n_bids-1)^(-1/alpha)*1/gamma_1p1oa*
             ((n_bids-1)*(gamma_1p1oa/mu*winning_bid)^alpha)^(1/alpha-1))
  }

  winning_bid + 1/alpha*(mu/gamma_1p1oa)*(n_bids-1)^(-1/alpha)*
    pgamma((n_bids-1)*(1/(mu/gamma_1p1oa)*winning_bid)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(n_bids-1)*(1/(mu/gamma_1p1oa)*winning_bid)^alpha)
}

vf__bid_function_fast = Vectorize(FUN = f__bid_function_fast,vectorize.args = "winning_bid")
