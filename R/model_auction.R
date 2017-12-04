#' Suite of functions to estimate private-value auction models
#'
#'
#' @param cost XXnumber of observations to draw
#' @param num_bids XXnon-negative alpha parameter of the beta distribution
#' @param mu XXnon-negative beta parameter of the beta distribution
#' @param alpha XXnon-negative beta parameter of the beta distribution
#' @param gamma_1p1oa XXnon-negative beta parameter of the beta distribution
#'
#' @details The Beta distribution with parameters \eqn{a} and \eqn{b} has
#' density:
#'
#' \deqn{
#'     \Gamma(a+b)/(\Gamma(a)\Gamma(b))x^(a-1)(1-x)^(b-1)
#' }
#'
#' for \eqn{a > 0}, \eqn{b > 0} and \eqn{0 \le x \le 1}.
#'
#' @examples
#' # Draw from beta distribution with parameters a = 1 and b = 3
#' beta_plot(a = 1, b = 3)
#'
#' @seealso \code{\link{rbeta}}, \code{\link{geom_density}}
#'
#'
#' @export
auction_v3 <- function(dat = NULL,
                       winning_bid = NULL, number_of_bids = NULL,
                       init_privatevalue_mu = NULL,
                       init_privatevalue_a = NULL,
                       init_control = NULL,
                       init_common_sd = NULL,
                       common_distributions = NULL,
                       num_cores = NULL
                       ) {
  # Initialize environment
  num_cores = auction_v3__init_env(num_cores=num_cores)

  # Validate distributions requested for unobserved heterogeneity
  common_distributions = auction_v3__check__common_distrib(common_distributions = common_distributions)

  # Validate input data
  auction_v3__check_input_data(dat = dat,
                               colName_price = winning_bid, colName_num = number_of_bids
                               )

  # Prepare initial guesses
  vecInitGuess = auction_v3__check_init_guess(dat = dat,
                                               colName_price = winning_bid, colName_num = number_of_bids,
                                               init_privatevalue_mu = init_privatevalue_mu,
                                               init_privatevalue_a = init_privatevalue_a,
                                               init_control = init_control,
                                               init_common_sd = init_common_sd
                                               )

  # Prepare control parameters for numerical solver
  conv_ctrl = auction_v3__get_conv_ctrl(vecInitGuess = vecInitGuess)

  # Set up parallelization of numerical solver
  cl = parallel::makeCluster(num_cores)
  # do we need to include ?
  #   do.call()
  #   dgamma()
  #   dlnorm()
  #   dweibull
  parallel::clusterExport(cl,
                varlist=c("vf__bid_function_fast__v4",
                          "vf__w_integrand_z_fast__v4",
                          "f__funk__v4",
                          "auction_v3__get_unobs_params",
                          "auction__get_distrib_params__gamma",
                          "auction__get_distrib_params__lognorm",
                          "auction__get_distrib_params__weibull")
                )
  # Run
  run_result = list()
  for (funcName in common_distributions) {
    sFuncName = as.character(funcName)

    # Run
    print(paste("Running |", sFuncName))
    #   Build function call parameter
    listFuncCall = list(funcName = sFuncName,
                        funcID = auction_v3__get_id_distrib(sFuncName = sFuncName),
                        argList = list())
    #   Run
    print("Fix dat_X=dat[['x_terms']] to a more generic solution")
    run_result[[sFuncName]] = stats::optim(par=vecInitGuess,
                                    fn=f__ll_parallel__v4,
                                    control=conv_ctrl,
                                    dat_price=dat[[winning_bid]],
                                    dat_num=dat[[number_of_bids]],
                                    dat_X=dat[['x_terms']],
                                    listFuncCall=listFuncCall,
                                    cl=cl)
    print(paste("End of run |", sFuncName))
  }
  # Release resources
  parallel::stopCluster(cl)
  # Prepare output
  res = auction_v3__output_org(run_result)
  return(res)
}

auction_v3__output_org <- function(run_result) {

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

    df['pv_weibull_mu', sFuncName] = run_result[[sFuncName]]$par[idxList$pv_weibull_mu]
    df['pv_weibull_a', sFuncName] = run_result[[sFuncName]]$par[idxList$pv_weibull_a]
    df['unobs_hetero__std_dev', sFuncName] = run_result[[sFuncName]]$par[idxList$unobs_dist_param]
    for (iX in 1:(1+length(run_result[[sFuncName]]$par)-idxList$x_terms__start)) {
      df[sprintf('X%d', iX), sFuncName] = run_result[[sFuncName]]$par[(idxList$x_terms__start+iX-1)]
    }

  }
  return(df)
}

auction__generate_data <- function(obs = 200) {
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

  return( list(
    price = y,
    num = n,
    x_terms = as.matrix( cbind( log(x1), log(x2) ) )
  ) )
}

auction__gen_err_msg <- function(res) {
  # Goal: Print out an error message and then stops execution of the main script

  errMsg = paste('\n\tError Code=', res['err_code'], '\n\tError msg=', res['err_msg'], sep='')
  stop(errMsg)
}

auction_v3__init_env <- function(num_cores) {
  # Goal:
  #   - Load all required packages, stop execution is any packages are missing
  #   - Check number of cores requested

  # Load required packages
  auction_v3__load_packages()
  # Check number of cores requested
  num_cores = auction_v3__check__num_cores(num_cores = num_cores)
  return(num_cores)
}

auction_v3__load_packages <- function () {
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

auction_v3__check__num_cores <- function(num_cores) {
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

auction_v3__check__common_distrib <- function(common_distributions) {
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

auction_v3__check_input_data <- function(dat, colName_price, colName_num) {
  # Goal: Make sure the input data has required columns

  if (mode(dat) == "list" && is.character(colName_price) && is.character(colName_num)) {
    # Get list of required column names
    colList_req = c(colName_price, colName_num)

    # Get list of column names
    if (is.data.frame(dat)) {
      colList = colnames(dat)
    } else {
      colList = names(dat)
    }

    # Make sure we aren't missing the required columns
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
    }
  } else {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'dat' and/or associated column names"
    auction__gen_err_msg(res)
  }
}

auction_v3__check_init_guess <- function(dat = dat,
                                         colName_price, colName_num,
                                         init_privatevalue_mu,
                                         init_privatevalue_a,
                                         init_control,
                                         init_common_sd) {

  # Find number of X parameters within dat
  #   Total number of columns - column price - column number of bids
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



  if (is.null(init_privatevalue_mu)) {
    x0[idxList$pv_weibull_mu] = def_pv_mu
  } else if (is.numeric(init_privatevalue_mu)) {
    x0[idxList$pv_weibull_mu] = init_privatevalue_mu
  } else {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'init_privatevalue_mu'"
    auction__gen_err_msg(res)
  }
  if (is.null(init_privatevalue_a)) {
    x0[idxList$pv_weibull_a] = def_pv_a
  } else if (is.numeric(init_privatevalue_a)) {
    x0[idxList$pv_weibull_a] = init_privatevalue_a
  } else {
    res = list()
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'init_privatevalue_a'"
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

  return(x0)
}

auction_v3__get_conv_ctrl <- function(vecInitGuess) {
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

auction_v3__get_id_distrib <- function(sFuncName) {
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

auction_v3__get_unobs_params <- function(distrib_std_dev, id_distrib) {
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

f__ll_parallel__v4 = function(x0, dat_price, dat_num, dat_X, listFuncCall, cl){
  # From iteration to iteration, only x0 is changing

  # update } else if ( sum ( ... ) ) {

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
    v__w = dat_price / v__h

    # Set E(X) = 1 for UnObserved distribution
    listFuncCall$argList = auction_v3__get_unobs_params(
      distrib_std_dev = x0[listIdx$unobs_dist_param],
      id_distrib = listFuncCall$funcID)

    # Run
    v__f_w = parallel::parApply(cl = cl,
                      X = cbind(v__w, dat_num, x0[listIdx$pv_weibull_mu], x0[listIdx$pv_weibull_a], v__gamma_1p1opa),
                      MARGIN = 1,
                      FUN = f__funk__v4,
                      listFuncCall=listFuncCall
    )

    # Return output
    ### print(paste("v__f_w", -sum(log(v__f_w / v__h))))
    ### print(func__prob_distrib)
    v__f_y = v__f_w / v__h
    return(-sum(log(v__f_y)))
  }
}

f__funk__v4 = function(data_vec, listFuncCall){
  val = stats::integrate(vf__w_integrand_z_fast__v4, w_bid=data_vec[1],
                  num_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], listFuncCall=listFuncCall,
                  lower=0, upper=Inf, abs.tol = 1e-10)
  if(val$message != "OK")
    stop("Integration failed.")
  return(val$value)
}

vf__w_integrand_z_fast__v4 = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, listFuncCall){

  # Get "x"
  b_z = vf__bid_function_fast__v4(price=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z

  # Add "x"
  listFuncCall$argList$x = u_z

  #Run
  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    do.call(
      match.fun(listFuncCall$funcName),
      listFuncCall$argList
    )
  ### dlnorm(u_z, meanlog=(-param_u^2*1/2), sdlog = param_u) # Note: can swap for different distributions

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}

f__bid_function_fast__v4 = function(price, num_bids, mu, alpha, gamma_1p1oa){

  if (exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*price)^alpha) == 0) {
    return(price + mu/alpha*(num_bids-1)^(-1/alpha)*1/gamma_1p1oa*
             ((num_bids-1)*(gamma_1p1oa/mu*price)^alpha)^(1/alpha-1))
  }

  price + 1/alpha*(mu/gamma_1p1oa)*(num_bids-1)^(-1/alpha)*
    stats::pgamma((num_bids-1)*(1/(mu/gamma_1p1oa)*price)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*price)^alpha)
  # Check gamma(1/alpha) part
}
vf__bid_function_fast__v4 = Vectorize(FUN = f__bid_function_fast__v4,vectorize.args = "price")


