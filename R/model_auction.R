#' Estimate private-value auction models
#'
#'
#' @param dat List containing the winning bids, number of bids, and \code{X} variables that describe the data.
#' @param winning_bid In list \code{dat}, the key whose value is a vector that holds the winning bids.
#' @param n_bids In list \code{dat}, the key whose value is a vector that holds the number of bids.
#' @param init_mu Value for \code{mu} for initial guess of the private value distribution.
#' @param init_alpha Value for \code{alpha} for initial guess of the private value distribution.
#' @param init_control Value(s) for the intial X terms used for the initla guess of the unobserved heterogeneity.
#' @param init_common_sd Value for \code{sigma} for initial guess of the private value distribution.
#' @param common_distributions Which distributions to test for modeling the unobserved heterogeneity.
#' @param num_cores The number of cores for running the model in parallel.
#'
#' @details This function attempts to describe auction data as a combination of private value data and unobserved heterogeneity.
#' Private values are modeled using the weibull distribution, while the unobserved heterogeneity can be modeled based on user input.
#' The function attempts to describe the data by iterative optimization of a minimization function, with the best fit reported
#' via the smallest standard deviation. X terms describing additional characterstics of the data can be supplied;
#' the probability density functions to use, and intial values for the shapes of the distrubtions, and intials values for the X
#' terms can be supplied, though they are optional. Initial, sample data can be supplied via the
#' \code{\link{auction_generate_data}} function. Otherwise, a minimum of 3 bids is required.
#'
#' \code{dat} must be a list having three keys, the names of which can be user-definted. The first key, \code{price}, is a numeric vector
#' describing the winniding bids. The second key, \code{num}, descibes the number of bids. And the last key, \code{x_terms}, is a numeric matrix whose
#' values describe the characteristics of the data. The key names must be supplied as the \code{winning_bid} and \code{n_bids} parameters.
#' VERIFY! The \code{x_terms} are determined by all remaining numeric columns in \code{dat}.
#'
#' Modeling for the unobserved heterogeneity is controled by the \code{coimmon_distributions}. This is eitehr a string or vector of stings that
#' state which distrubtions to use for modeling: \code{dlnorm}, \code{weibull}, and \code{weibull}. This parameter is optional; the log normal
#' distribution \code{dlnrom} will be used if not supplied.
#'
#' \code{init_mu}, \code{init_alpha}, \code{init_common_sd}, and \code{init_control} are all optional parameters. If not supplied, XXX.
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
#'
#'
#' @return For each of the distributions speicifed in \code{common_distributions}, ...
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
model_auction <- function(dat = NULL,
                       winning_bid = NULL, n_bids = NULL,
                       init_mu = NULL,
                       init_alpha = NULL,
                       init_control = NULL,
                       init_common_sd = NULL,
                       init_params = NULL,
                       common_distributions = NULL,
                       num_cores = 1
                       ) {
  # Initialize environment
  num_cores = auction__init_env(num_cores=num_cores)

  # Validate distributions requested for unobserved heterogeneity
  common_distributions = auction__check__common_distrib(common_distributions = common_distributions)

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
                                           init_control = init_control,
                                           init_common_sd = init_common_sd,
                                           init_params=init_params
                                           )

  # Prepare control parameters for numerical solver
  conv_ctrl = auction__get_conv_ctrl(vecInitGuess = vecInitGuess)

  # Set up parallelization of numerical solver
  cl = parallel::makeCluster(num_cores)
  # do we need to include ?
  #   do.call()
  #   dgamma()
  #   dlnorm()
  #   dweibull
  parallel::clusterExport(cl,
                varlist=c("vf__bid_function_fast",
                          "vf__w_integrand_z_fast",
                          "f__funk",
                          "auction__get_unobs_params",
                          "auction__get_distrib_params__gamma",
                          "auction__get_distrib_params__lognorm",
                          "auction__get_distrib_params__weibull"),
                envir = environment(model_auction)
                )

  # Run
  run_result = list()
  for (funcName in common_distributions) {
    sFuncName = as.character(funcName)

    # Run
    print(paste("Running |", sFuncName))
    #   Build function call parameter
    listFuncCall = list(funcName = sFuncName,
                        funcID = auction__get_id_distrib(sFuncName = sFuncName),
                        argList = list())
    #   Run
    run_result[[sFuncName]] = stats::optim(par=vecInitGuess,
                                           fn=f__ll_parallel,
                                           control=conv_ctrl,
                                           dat__winning_bid=dat[[winning_bid]],
                                           dat__n_bids=dat[[n_bids]],
                                           dat_X=dat[ ! names(dat) %in% c(winning_bid, n_bids) ],
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

#' Generate example data for running \code{\link{model_auction}}
#'
#'
#' @param obs Number of observations to draw
#'
#' @details This function generates example data for feeding into model_auction(). Specifically, the
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
#' @seealso \code{\link{model_auction}}
#'
#'
#' @export
auction_generate_data <- function(obs = 200) {
  # For testing purposes, we will generate sample data

  # ensure that obs is integer greater than 0
  set.seed(301)
  # data = # Generate some data
  # y, n, x1, x2: positive
  # n: discrete and > 1
  # y is some function of n, x1, x2

  w = stats::rlnorm(obs)
  x1 = stats::rlnorm(obs) + 0.5*w
  x2 = 0.1*stats::rlnorm(obs) + 0.3*w
  e = 2*stats::rlnorm(obs)
  n = sample(2:10, obs, replace=TRUE)
  y = 10 - 0.5*n + x1 + x2 + e

  # return( list(
  #   price = y,
  #   num = n,
  #   x_terms = as.matrix( cbind( log(x1), log(x2) ) )
  # ) )
  return( data.frame(cbind(y, n, x1, x2)) )
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
auction__check__common_distrib <- function(common_distributions) {
  if (is.null(common_distributions)) {
    common_distributions = 'dlnorm'
  } else if (is.character(common_distributions)) {
    for (distrib in common_distributions) {
      if ((distrib != 'dlnorm') && (distrib != 'dweibull') && (distrib != 'dgamma')) {
        res = list()
        res['err_code'] = 2
        res['err_msg'] = paste("Invalid input for 'common_distributions' | Entry '", distrib, "' is invalid", sep='')
        auction__gen_err_msg(res)
      }
    }
  } else {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'common_distributions' | Must be string or vector of strings"
    auction__gen_err_msg(res)
  }
  return(common_distributions)
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
                                      init_control,
                                      init_common_sd,
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

    if (is.null(init_common_sd)) {
      x0[idxList$unobs_dist_param] = def_pv_mu
    } else if (is.numeric(init_common_sd)) {
      x0[idxList$unobs_dist_param] = init_common_sd
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_common_sd'"
      auction__gen_err_msg(res)
    }

    if (is.null(init_control)) {
      x0[idxList$x_terms__start:length(x0)] = def_x
    } else if (is.numeric(init_control)) {
      x0[idxList$x_terms__start:length(x0)] = init_control
    } else {
      res = list()
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'init_control'"
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
  # # Given std dev and E(X) = 1, calculate meanlog and sdlog
  # tmp = log(1+distrib_std_dev^2)
  # return(list(sdlog=sqrt(tmp), meanlog=-1/2*tmp))

  return(list( meanlog=(-distrib_std_dev^2*1/2), sdlog = distrib_std_dev ))
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
  # From iteration to iteration, only x0 is changing

  # Get position indices
  listIdx = auction__x0_indices()

  h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[2]-1)]
  v__h = exp( colSums( h * t(dat_X) ) )

  if (x0[listIdx$unobs_dist_param] <= 0.1) {
    return(-Inf) # Check that these hold at estimated values
  } else if ( sum (x0[listIdx$pv_weibull_mu] <= 0 ) > 0 ) {
    return(-Inf)
  } else if ( sum( x0[listIdx$pv_weibull_a] <= 0.01 ) > 0) {
    return(-Inf)
  } else {
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/x0[listIdx$pv_weibull_a])
    v__w = dat__winning_bid / v__h

    # Set E(X) = 1 for UnObserved distribution
    listFuncCall$argList = auction__get_unobs_params(
      distrib_std_dev = x0[listIdx$unobs_dist_param],
      id_distrib = listFuncCall$funcID)
    # Run
    v__f_w = parallel::parApply(cl = cl,
                                X = cbind(v__w,
                                          dat__n_bids,
                                          x0[listIdx$pv_weibull_mu],
                                          x0[listIdx$pv_weibull_a],
                                          v__gamma_1p1opa),
                                MARGIN = 1,
                                FUN = f__funk,
                                listFuncCall=listFuncCall
    )
    # Return output
    v__f_y = v__f_w / v__h
    return(-sum(log(v__f_y)))
  }
}
f__funk = function(data_vec, listFuncCall){
  val = stats::integrate(vf__w_integrand_z_fast, w_bid=data_vec[1],
                  n_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], listFuncCall=listFuncCall,
                  lower=0, upper=Inf, abs.tol = 1e-10)
  if(val$message != "OK")
    stop("Integration failed.")
  return(val$value)
}
vf__w_integrand_z_fast = function(z, w_bid, n_bids, mu, alpha, gamma_1p1oa,
                                  listFuncCall){

  # Get "x"
  b__winning_bid = vf__bid_function_fast(winning_bid=z, n_bids=n_bids,
                              mu=mu, alpha=alpha, gamma_1p1oa)
  u__winning_bid = w_bid/b__winning_bid

  # Add "x"
  listFuncCall$argList$x = u__winning_bid

  #Run
  vals = n_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-n_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b__winning_bid*
    do.call(
      match.fun(listFuncCall$funcName),
      listFuncCall$argList
    )
  ### dlnorm(u_z, meanlog=(-param_u^2*1/2), sdlog = param_u) # Note: can swap for different distributions

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
    stats::pgamma((n_bids-1)*(1/(mu/gamma_1p1oa)*winning_bid)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(n_bids-1)*(1/(mu/gamma_1p1oa)*winning_bid)^alpha)
  # Check gamma(1/alpha) part
}
vf__bid_function_fast = Vectorize(FUN = f__bid_function_fast,vectorize.args = "winning_bid")








#######################################################
# Load Data
#######################################################
set.seed(301)
# data = # Generate some data
# y, n, x1, x2: positive
# n: discrete and > 1
# y is some function of n, x1, x2

obs = 200
num_cores = 3

w = rlnorm(obs)
x1 = rlnorm(obs) + .5*w
x2 = .1*rlnorm(obs) + .3*w
e = 2*rlnorm(obs)
n = sample(2:10, obs, replace=TRUE)
y = 10 - .5*n + x1 + x2 + e
data = data.frame(cbind(y, n, x1, x2))
plot(n, y)

v.y = data$y
v.n = data$n
# m.h_x = as.matrix(cbind(log(data$x1),log(data$x2)))
m.h_x = as.matrix(cbind(data$x1,data$x2))

# inital parameter guess
x0 =  c(8, 2, .5, .4, .6)


# My estimation routine
# $par
# [1] 15.4977557  3.8299721  0.2145705  0.1577268  0.1068006
# $value
# [1] 508.4587
# $counts
# function gradient
# 495       NA


sdlog = x0[3]
meanlog = -1/2*sdlog^2
initial_guess_value =  sqrt((exp( sdlog^2 ) - 1) * exp( 2*meanlog + sdlog^2 ))

initial_guess_value = sdlog




init_params = c(x0[1], x0[2], initial_guess_value, x0[-c(1,2,3)])




data = data.frame("price" = v.y, "num" = v.n, "x1_name" = m.h_x[,1], "x2_name" = m.h_x[,2])
for (num_cores in c(2)) {
  t2 = system.time({
  res = model_auction(dat = data, winning_bid = 'price',
                n_bids = 'num',
                init_mu = x0[1],
                init_alpha = x0[2],
                init_common_sd = initial_guess_value,
                init_control = x0[-c(1,2,3)],
                init_params = init_params,
                num_cores = num_cores)
  })[[3]]
  print(paste("obs=", obs))
  print(paste("num_core=", num_cores))
  print(paste("runtime=", t2))

  print(res)
}
