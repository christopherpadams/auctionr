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
auction_v2 <- function(dat = NULL,
                       winning_bid = NULL,
                       number_of_bids = NULL,
                       initial_guess = NULL,
                       num_cores = NULL,
                       func_list__unobs_distrib = NULL
                       ) {
  #
  #
  #   winning_bid
  #     column name within dat that represents price or winning bid
  #
  #   number_of_bids
  #     column name within dat that represents number of bids


  # Initialize environment
  #   Initialize
  res = auction__init_env(num_cores=num_cores, func_list__unobs_distrib=func_list__unobs_distrib)
  #   Unpack 'res'
  num_cores = res$num_cores
  func_list__unobs_distrib = res$func_list__unobs_distrib
  res = res$res

  # Validate that data, winning_bid and number_of_bids all match
  res = auction__input_data_consistency(res=res,
                                        dat=dat,
                                        colName_price=winning_bid,
                                        colName_num=number_of_bids)

  #   Part 2) Gather components
  #     we shouldn't need to do this here
  #       when we use optim, we want to change input from
  #           optim( ... , dat_price=dat_price, dat_num=dat_num, ... )
  #       to
  #           optim( ... , dat_price=dat[[winning_bid]], ... )
  #v__y = dat$v__y
  #v__n = dat$v__n
  #m__h_x = dat$m__h_x
  dat_price=dat[[winning_bid]]
  dat_num=dat[[number_of_bids]]
  dat_X=dat[['x_terms']]


  # Prepare initial guess
  x0 = auction__conv__init_pos(res=res,
                               dat=dat,
                               x0=initial_guess)

  # Set up parallelization
  cl = makeCluster(num_cores)
  clusterExport(cl,varlist=c("vf__bid_function_fast__v3",
                             "vf__w_integrand_z_fast__v3",
                             "f__funk__v3"))
  #   Get convergence control parameters
  conv_ctrl = auction__conv__control_params__set(res=res,
                                                 x0=x0)

  #   Ensure initial guess, x0, is a vector not a list
  x0 = auction__conv__initial_guess__set(res=res,
                                         x0=x0)

  # Run
  run_result = list()
  for (funcName in names(func_list__unobs_distrib)) {
    sFuncName = as.character(funcName)

    # Initialize results
    run_result[[sFuncName]] = NULL
    # Prepare inputs
    func__unobs_distrib = list(funcName=sFuncName, argList=func_list__unobs_distrib[[sFuncName]])

    # Run
    print(paste("Running |", sFuncName))

    run_result[[sFuncName]] = optim(par=x0,
                                    fn=f__ll_parallel__v3,
                                    control=conv_ctrl,
                                    dat_price=dat_price,
                                    dat_num=dat_num,
                                    dat_X=dat_X,
                                    func__prob_distrib=func__unobs_distrib,
                                    cl=cl)

    print(paste("End of run |", sFuncName))
  }


  # Release resources
  stopCluster(cl)
  # Inspect result

  print(run_result)

  # For each distribution function representing unobserved heterogeneity
  #   Find parameter value
  #     Get index of x0 vector
  idxList = auction__x0_indices()

  run_result2 = list()
  for (funcName in names(run_result)) {
    sFuncName = as.character(funcName)

    # Get parameter that was being changed
    paramName = auction__unobs_dist__paramName(
      func__prob_distrib=list(funcName=sFuncName, argList=func_list__unobs_distrib[[sFuncName]])
      )
    # Get final parameter value
    paramVal = run_result[[sFuncName]]$par[idxList$unobs_dist_param]

    # Use auction__unobs_dist__exp_val_1() to get other parameter plus its value
    #   Build input, func__prob_distrib
    func__prob_distrib = list(funcName=sFuncName, argList=list())
    func__prob_distrib$argList[[paramName]] = paramVal
    #   Get output
    func__prob_distrib = auction__unobs_dist__exp_val_1(func__prob_distrib, paramVal)

    # Get standard deviation of unobserved heterogeneity
    if (sFuncName == 'dgamma') {

    } else if (sFuncName == 'dlnorm') {
      unobs_std = auction__unobs_dist__std__lognorm(
        param_meanlog=func__prob_distrib$argList$meanlog,
        param_sdlog=func__prob_distrib$argList$sdlog
        )
    } else if (sFuncName == 'dweibull') {

    }

    run_result2[[sFuncName]] = list()
    run_result2[[sFuncName]]$invLogLikelyhood = run_result[[sFuncName]]$value
    run_result2[[sFuncName]]$argList = func__prob_distrib$argList
    run_result2[[sFuncName]]$std = unobs_std
  }
  # print(run_result2)


  # If inspection is fine, add to res
  res$result = run_result2
  # Return result
  return(res)
}


auction__generate_data <- function(obs = 200) {
  # For testing purposes, we will generate sample data

  # ensure that obs is integer greater than 0
  set.seed(301)
  # data = # Generate some data
  # y, n, x1, x2: positive
  # n: discrete and > 1
  # y is some function of n, x1, x2

  w = rlnorm(obs)
  x1 = rlnorm(obs) + 0.5*w
  x2 = 0.1*rlnorm(obs) + 0.3*w
  e = 2*rlnorm(obs)
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
auction__init_env <- function (num_cores, func_list__unobs_distrib) {
  # Goal:
  #   - Prepare results/error handling
  #   - Load all required packages, stop execution is any packages are missing
  #   - Check number of cores requested
  #   - Check parameters related to the unobserved heterogeneity

  # Initialize output and error code + message
  res = list(result = -1, err_code = 0, err_msg = "")
  # Load required libraries
  res = auction__load_packages_v2(res)
  # Check the number of cores requested
  #   result = list(res=res, num_cores=num_cores)
  res = auction__check_input__num_cores__v2(res, num_cores)
  #     unpack
  num_cores = res$num_cores
  res = res$res
  # Check PDFs specified
  res = auction__check_input__func_list__unobs_distrib(func_list__unobs_distrib, res)
  #   unpack
  func_list__unobs_distrib = res$func_list
  res = res$res

  return( list(
    res = res,
    num_cores = num_cores,
    func_list__unobs_distrib = func_list__unobs_distrib
  ) )
}
auction__load_packages_v2 <- function (res) {
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
    res['err_code'] = 1
    res['err_msg'] = paste0("Unable to load the following packages: ",
                            paste(listMissingPackages, collapse=','))

    auction__gen_err_msg(res)
  }
  return(res)
}
auction__check_input__num_cores__v2 <- function(res, num_cores) {
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
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'num_cores' | Must be natural number"
    auction__gen_err_msg(res)
  } else {
    # Check cores available
    num_cores__avail = detectCores(res)
    if ( num_cores > num_cores__avail ) {
      print(paste0("Warning: You have requested ", num_cores,
                   " cores but only have ", num_cores__avail,
                   " cores available"))

      # Adjust number of workers we will request
      num_cores = num_cores__avail
      print(paste0("\tSetting # parallel workers to", num_cores))
    }
  }
  return(list(res=res, num_cores=num_cores))
}
auction__input_data_consistency <- function(res, dat, colName_price, colName_num) {
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
      res['err_code'] = 2
      res['err_msg'] = paste0("Unable to find the following columns within input data: ",
                              paste(listMissingColName, collapse=','))
      auction__gen_err_msg(res)
    }
  } else {
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'dat' and/or associated column names"
  }
  if (res['err_code'] != 0) {
    auction__gen_err_msg(res)
  }
  return(res)
}

auction__conv__init_pos <- function(res, dat, x0) {
  # Goal: Validate structure or entries within x0
  #           x0 is the initial guess or initial position for
  #           the numerical solver
  #             -> optim()

  # Check if user-supplied 'x0' is valid
  #     Length of x0 must be 2+1+ # of Xs
  #     Width of dat must be 2 + # of Xs
  #
  #   First need to know how many X parameters are within data
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
  nX_dat = nParams_dat - 2

  #   Now check if x0 is valid
  if (! is.null(x0)) {
    if (is.numeric(x0)) {
      nParams_x0 = length(x0)
    } else if (is.list(x0)) {
      nParams_x0 = 0
      for (colName in names(x0)) {
        nCol = dim(x0[[colName]])[2]
        if (is.null(nCol)) {
          nCol = 1
        }
        nParams_x0 = nParams_x0 + nCol
      }
    } else {
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'initial_guess' | Must be numeric vector or list"
      auction__gen_err_msg(res)
    }

    nX_x0 = nParams_x0 - 3
    if (nX_x0 < nX_dat) {
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'initial_guess' | Too few parameters given"
      auction__gen_err_msg(res)
    }
  } else {
    # No initial guess given, generate
    x0 = auction__conv__init_pos__gen(dat__nParams=nParams_dat)
  }

  return(x0)
}
auction__conv__init_pos__gen <- function(dat__nParams) {
  # Goal: Generate x0
  #           x0 is the initial guess or initial position for
  #           the numerical solver
  #             -> optim()

  # Currently, fix initial guess
  pv_weibull_mu = 8
  pv_weibull_a = 2

  unobs_dist_param = 0.5

  # dat = "price" | "num of bids" | "x1" | "x2" | ...
  #   # X columns = # Columns - 2
  h_x = rep(0.5, (dat__nParams-2) )

  return( list(
    pv_weibull_mu = pv_weibull_mu,
    pv_weibull_a = pv_weibull_a,
    unobs_dist_param = unobs_dist_param,
    x_terms = h_x
    ) )
}

# auction__conv__control_params__set <- function(res, x0) {
#   # Define convergence criteria
#   #  Max iterations
#   #  Step-size of parameters that may be adjusted
#
#   max_iterations = 2000
#
#
#
#
#   val__pv_weibull_mu = NULL
#   idx__pv_weibull_mu = 1
#   val__pv_weibull_a = NULL
#   idx__pv_weibull_a = 2
#   val__unobs_dist_param =NULL
#   idx__unobs_dist_param =3
#   val__x_terms = NULL
#   idx__x_terms__start = 4
#
#
#
#
#
#   pv_weibull_mu__initial_guess = 1
#   pv_weibull_a__initial_guess = 0.1
#   unobs_dist_param__initial_guess = 1
#   h_x__1__initial_guess = 0.1
#   h_x__i__initial_guess = 1
#
#   if (is.numeric(x0)) {
#     # x0 should be of the form
#     #   c( [pv_weibull_mu] , [pv_weibull_a] , [unobs_dist_param] , [h_x 1] , [h_x 2] , ... )
#     return( list(
#       maxit = max_iterations,
#       parscale = c(
#         pv_weibull_mu__initial_guess,
#         pv_weibull_a__initial_guess,
#         unobs_dist_param__initial_guess,
#         h_x__1__initial_guess,
#         rep(h_x__i__initial_guess, length(x0) - 4)
#         )
#     ) )
#   } else if (is.list(x0)) {
#
#     x0_step = c()
#     step_default = 1
#
#     for ( x0_key in names(x0) ) {
#
#       if (length(grep('unobs', x0_key))>=1) {
#         x0_step = c(x0_step, unobs_dist_param__initial_guess)
#       } else if (length(grep('weibull', x0_key))>=1) {
#         if (length(grep('_mu', x0_key))==1) {
#           x0_step = c(x0_step, pv_weibull_mu__initial_guess)
#         } else if (length(grep('_a', x0_key))>=1) {
#           x0_step = c(x0_step, pv_weibull_a__initial_guess)
#         } else {
#           x0_step = c(x0_step, step_default)
#         }
#       } else if (length(grep('x_term', x0_key))>=1) {
#         nCol = length(x0[[x0_key]])
#         if (nCol == 1) {
#           x0_step = c(x0_step, h_x__1__initial_guess)
#         } else {
#           x0_step = c(x0_step, h_x__1__initial_guess, rep(h_x__i__initial_guess, nCol-1))
#         }
#
#       } else {
#         x0_step = c(x0_step, step_default)
#       }
#
#     }
#
#     return( list(
#       maxit = max_iterations,
#       parscale = x0_step
#       ) )
#
#   } else {
#     res['err_code'] = 2
#     res['err_msg'] = "Unexpected error, x0 is niether numeric vector nor list"
#     auction__gen_err_msg(res)
#   }
# }
auction__conv__control_params__set <- function(res, x0) {
  # Define convergence criteria
  #  Max iterations
  #  Step-size of parameters that may be adjusted

  # define maxit
  max_iterations = 2000

  # ready to define step sizes
  #   x0 step size list should be of the form
  #       c( [pv_weibull_mu] , [pv_weibull_a] , [unobs_dist_param] , [h_x 1] , [h_x 2] , ... )
  #   get indices and step sizes
  listIdx = auction__x0_indices()
  listStep = auction__x0_stepsizes()
  #   get number of parameters in input x0 and number of X terms
  if (is.numeric(x0)) {
    # get number of parameters in input x0
    nParam = length(x0)
    # get number of X terms
    #   X terms are always in the last positions of the numeric vector x0
    nParam_X = length(x0) - listIdx$x_terms__start + 1
  } else if (is.list(x0)) {
    nParam = 0
    nParam_X = 0
    for ( x0_key in names(x0) ) {
      nParam = nParam + length(x0[[x0_key]])
      if (length(grep('x_term', x0_key))>=1) {
        nParam_X = nParam_X + length(x0[[x0_key]])
      }
    }
    if (nParam_X == 0) {
      res['err_code'] = 2
      res['err_msg'] = "Unexpected error, x0 is list but no X terms found"
      auction__gen_err_msg(res)
    }
  } else {
    res['err_code'] = 2
    res['err_msg'] = "Unexpected error, x0 is niether numeric vector nor list"
    auction__gen_err_msg(res)
  }
  # build vector containg step sizes for just X terms
  x0_steps__x_terms = c(listStep$x_terms__start, rep(listStep$x_terms__other, nParam_X-1))
  # initialize vector to contain all step sizes
  x0_steps = rep(-1, nParam)
  # fill vector with all step sizes
  if (is.numeric(x0)) {
    # We assume x0 holds the correct parameters and is in the correct ordering,
    #     c( [pv_weibull_mu] , [pv_weibull_a] , [unobs_dist_param] , [h_x 1] , [h_x 2] , ... )
    x0_steps[listIdx$pv_weibull_mu] = listStep$pv_weibull_mu
    x0_steps[listIdx$pv_weibull_a] = listStep$pv_weibull_a

    x0_steps[listIdx$unobs_dist_param] = listStep$unobs_dist_param

    x0_steps[listIdx$x_terms__start:(listIdx$x_terms__start + nParam_X - 1)] = x0_steps__x_terms

  } else if (is.list(x0)) {

    for ( x0_key in names(x0) ) {

      if (length(grep('unobs', x0_key))>=1) {
        x0_steps[listIdx$unobs_dist_param] = listStep$unobs_dist_param
      } else if (length(grep('weibull', x0_key))>=1) {
        if (length(grep('_mu', x0_key))==1) {
          x0_steps[listIdx$pv_weibull_mu] = listStep$pv_weibull_mu
        } else if (length(grep('_a', x0_key))>=1) {
          x0_steps[listIdx$pv_weibull_a] = listStep$pv_weibull_a
        }
      } else if (length(grep('x_term', x0_key))>=1) {
        x0_steps[listIdx$x_terms__start:(listIdx$x_terms__start + nParam_X - 1)] = x0_steps__x_terms
      }

    }

  }

  if (any(x0_steps == -1)) {
    res['err_code'] = 2
    res['err_msg'] = "Unexpected error, x0_step has at least one -1 entry"
    if (is.numeric(x0)) {
      res['err_msg'] = paste(res['err_msg'], '| x0 is numeric')
    } else {
      res['err_msg'] = paste(res['err_msg'], '| x0 is list')
    }
    auction__gen_err_msg(res)
  }
  return( list(
    maxit = max_iterations,
    parscale = x0_steps
  ) )
}






auction__conv__initial_guess__set <- function(res, x0) {
  # x0 must be of the form
  #   c( [pv_weibull_mu] , [pv_weibull_a] , [unobs_dist_param] , [h_x 1] , [h_x 2] , ... )
  if (is.numeric(x0)) {
    # Assumed correct
    return(x0)
  } else if (is.list(x0)) {
    # Must convert to vector
    #   Define indices
    listIdx = auction__x0_indices()
    #   Initialize values to NULL
    val__pv_weibull_mu = NULL
    # idx__pv_weibull_mu = 1
    val__pv_weibull_a = NULL
    # idx__pv_weibull_a = 2
    val__unobs_dist_param =NULL
    # idx__unobs_dist_param =3
    val__x_terms = NULL
    # idx__x_terms__start = 4

    #   Initialize vector
    #     Find length
    nParam = 0
    for ( x0_key in names(x0) ) {
      nParam = nParam + length(x0[[x0_key]])
    }
    #     Initialize to -1
    x0_vec = rep(-1, nParam)
    #   Fill vector
    for ( x0_key in names(x0) ) {

      if (length(grep('unobs', x0_key))>=1) {

        if (is.null(val__unobs_dist_param)) {
          val__unobs_dist_param = x0[[x0_key]]
          x0_vec[listIdx$unobs_dist_param] = x0[[x0_key]]
        } else {
          res['err_code'] = 2
          res['err_msg'] = "Unexpected error, x0 contains multiple parameters for unobserved distrib function"
        }

      } else if ((length(grep('weibull', x0_key))>=1) && (length(grep('_mu', x0_key))==1)) {

          if (is.null(val__pv_weibull_mu)) {
            val__pv_weibull_mu = x0[[x0_key]]
            x0_vec[listIdx$pv_weibull_mu] = x0[[x0_key]]
          } else {
            res['err_code'] = 2
            res['err_msg'] = "Unexpected error, x0 contains multiple 'mu' parameters for private-value distrib function"
          }

      } else if ((length(grep('weibull', x0_key))>=1) && (length(grep('_a', x0_key))==1)) {

        if (is.null(val__pv_weibull_a)) {
          val__pv_weibull_a = x0[[x0_key]]
          x0_vec[listIdx$pv_weibull_a] = x0[[x0_key]]
        } else {
          res['err_code'] = 2
          res['err_msg'] = "Unexpected error, x0 contains multiple 'a' parameters for private-value distrib function"
        }

      } else if (length(grep('x_term', x0_key))>=1) {

        if (is.null(val__x_terms)) {
          val__x_terms = x0[[x0_key]]

          nParam_X = length(x0[[x0_key]])
          x0_vec[listIdx$x_terms__start:(listIdx$x_terms__start+nParam_X-1)] = x0[[x0_key]]

        } else {
          res['err_code'] = 2
          res['err_msg'] = "Unexpected error, x0 contains multiple 'X' parameters"
        }

      } else {
        res['err_code'] = 2
        res['err_msg'] = "Unexpected error, x0 contains unknown parameters"
      }

    }

    if (res['err_code'] != 0) {
      auction__gen_err_msg(res)
    } else if (
      is.null(val__pv_weibull_mu) || is.null(val__pv_weibull_a) || is.null(val__unobs_dist_param) || is.null(val__x_terms)
      ) {
      res['err_code'] = 2
      res['err_msg'] = "Unexpected error, x0 is list() and missing parameters"
      auction__gen_err_msg(res)
    } else {
      return(x0_vec)
    }

  } else {
    res['err_code'] = 2
    res['err_msg'] = "Unexpected error, x0 is niether numeric vector nor list"
    auction__gen_err_msg(res)
  }
}


auction__check_input__func_list__unobs_distrib <- function(func_list, res) {
  # Goal: Check parameters related to the unobserved heterogeneity

  sInputVar = "'func_list__unobs_distrib'"

  if (is.null(func_list)) {
    func_list = list(dlnorm=NULL)
  }

  if (! is.list(func_list) || length(func_list) == 0) {
    res['err_code'] = 2
    res['err_msg'] = paste("Invalid input for", sInputVar, '| Must be NULL or list() as per documentation')
    # res['err_msg'] = "Invalid input for 'pdf_list' | Must be NULL for default PDF or list as per documentation"

  } else {
    func_list2 = func_list

    # Get list of valid PDFs
    func_list__valid = auction__valid_opt__func_list()

    # Remove PDFs that are not valid
    func_list2 = func_list2[names(func_list2) %in% names(func_list__valid)]

    if (length(func_list2)==0) {
      res['err_code'] = 2
      res['err_msg'] = paste("Invalid probability distribution functions selected for", sInputVar)
      # res['err_msg'] = "Invalid PDFs specified for 'pdf_list'"
    } else {

      # Check to see if parameter requirements are met
      for (funcName in names(func_list2)) {
        I_valid = TRUE

        curFunc = func_list2[[as.character(funcName)]]

        # Remove parameters that are not valid
        if (is.list(curFunc)) {
          # Remove by parameter name
          curFunc = curFunc[names(curFunc) %in% names(func_list__valid[[as.character(funcName)]])]

          if (length(names(curFunc)) != 0){
            # Remove by parameter datatype
            for (paramName in names(curFunc)) {
              # Get parameter value being tested
              paramVal = curFunc[[as.character(paramName)]]
              #   If NULL, then remove
              #   Otherwise, if not boolean nor numeric, then wrap in quotes
              if (is.null(paramVal)) {
                # Remove and continue
                curFunc = curFunc[names(curFunc) != paramName]
                next
              } else if ((! is.logical(paramVal)) && (! is.numeric(paramVal))) {
                # Wrap in quotes
                paramVal = paste("'", paramVal, "'", sep="")
              }
              # Get required parameter data type
              paramReqType = func_list__valid[[as.character(funcName)]][[as.character(paramName)]][['type']]

              # Check
              checkStr = paste("is.", paramReqType, "(", paramVal, ")", sep="")
              I_check = eval(parse(text=checkStr))

              if (! I_check) {
                # Remove
                curFunc = curFunc[names(curFunc) != paramName]
              }

            }
          }
        }

        # Check for required parameters
        #   curPDF or length(curPDF)==0 is okay if no required parameters
        if (is.null(curFunc) || (is.list(curFunc))) {
          # Look for requirements
          req_list = c()
          validFunc = func_list__valid[[as.character(funcName)]]
          if (is.list(validFunc) && length(validFunc) > 0) {
            for (paramName in names(validFunc)) {

              if (validFunc[[as.character(paramName)]][['req']]) {
                req_list = c(req_list, paramName)
              }

            }
          }

          # Check to see if requirements are met
          if (length(req_list) > 0) {
            if (is.null(curFunc) || (length(names(curFunc))==0)) {
              I_valid = FALSE
            } else {

              I_req = reqList %in% names(curFunc)
              if (! all(I_req)) {
                I_valid = FALSE
                print(paste(
                  "Ignoring probability distribution function '", funcName, "' because missing required parameters: ",
                  paste(reqList[! I_req], collapse=','), sep = ""
                ))
              }

            }
          }

        } else {
          I_valid = FALSE
        }

        if (! I_valid) {
          func_list2 = func_list2[names(func_list2) != funcName]
        } else {
          if (length(names(curFunc))==0) {
            curFunc = list(NULL)
          }

          func_list2[[as.character(funcName)]] = curFunc

        }
      }

      if (length(func_list2)==0) {
        res['err_code'] = 2
        res['err_msg'] = paste('Invalid input for', sInputVar, '| No valid distribution function or missing required parameter(s)')
        # res['err_msg'] = paste("No valid probability distribution functions within", sInputVar, 'or missing required parameters')
        # res['err_msg'] = "No PDFs specified in 'pdf_list' have valid required parameters"
      }

    }

    func_list = func_list2
  }

  if (res['err_code'] != 0) {
    auction__gen_err_msg(res)
  }
  return(list(res=res,func_list=func_list))
}
auction__valid_opt__func_list <- function() {
  # Goal: Provide list of functions and associated parameters that are accepted for representing
  #       unobserved heterogenity
  #         functions = valid probability distribution fuction, all are 2-parameter functions
  #                     each function could  accept the argument 'log'

  # Parameters have two values, (1) whether they are required and (2) what datatype must they be
  #   If a parameter is not required, than (1) is NULL
  return(list(
    dgamma=list(
      shape=list(type='numeric',req=FALSE),
      rate=list(type='numeric',req=FALSE)
    ),
    dlnorm=list(
      meanlog=list(type='numeric',req=FALSE),
      sdlog=list(type='numeric',req=FALSE)
    ),
    dweibull=list(
      shape=list(type='numeric',req=TRUE),
      scale=list(type='numeric',req=FALSE)
    )
  ))
}



auction__unobs_dist__paramName <- function(func__prob_distrib) {
  if (is.null( names(func__prob_distrib$argList) )) {
    if (func__prob_distrib$funcName == 'dgamma') {
      paramName = 'rate'
    } else if (func__prob_distrib$funcName == 'dlnorm') {
      paramName = 'sdlog'
    } else if (func__prob_distrib$funcName == 'dweibull') {
      paramName = 'shape'
    }
  } else {
    paramName = names(func__prob_distrib$argList)[1]
  }
  return(paramName)
}

auction__unobs_dist__exp_val_1 <- function(func__prob_distrib, paramVal) {
  # Goal: Update parameters within 'func__prob_distrib' to ensure
  #       E(X)=1 for the probability distribution representing
  #       unobserved heterogeneity

  # Get name of parameter that we are focusing on
  paramName = auction__unobs_dist__paramName(func__prob_distrib=func__prob_distrib)

  # Check if we are to hold constant the value of the parameter we are focusing on,
  #   Afterwards, calculate value for other parameter to get E(x) = 1
  if (! is.null( names(func__prob_distrib$argList) )) {
    # Is the value being held constant?
    if (!is.null(func__prob_distrib$argList[[paramName]])) {
      paramVal = func__prob_distrib$argList[[paramName]]
    }
  }
  # Ready to calculate value for other parameter
  if (func__prob_distrib$funcName == 'dgamma') {

  } else if (func__prob_distrib$funcName == 'dlnorm') {
    # Calculate
    if (paramName == 'meanlog') {
      list_paramVal = auction__unobs_dist__exp_val_1__lognorm(param_meanlog=paramVal)
    } else {
      list_paramVal = auction__unobs_dist__exp_val_1__lognorm(param_sdlog=paramVal)
    }
    # Store values with $argList
    if (is.null( names(func__prob_distrib$argList) )) {
      func__prob_distrib$argList = list(
        meanlog=list_paramVal$param_meanlog,
        sdlog=list_paramVal$param_sdlog
      )
    } else {
      func__prob_distrib$argList$meanlog = list_paramVal$param_meanlog
      func__prob_distrib$argList$sdlog = list_paramVal$param_sdlog
    }

  } else if (func__prob_distrib$funcName == 'dweibull') {

  }
  return(func__prob_distrib)
}

auction__unobs_dist__exp_val_1__lognorm <- function(param_meanlog = NULL, param_sdlog = NULL) {
  # Given either meanlog or sdlog, set the other such that E(x) = 1
  #   E(x) = exp(mu) * exp(1/2 * s ^ 2)

  if (is.null(param_meanlog)) {
    param_meanlog = -1/2 * param_sdlog ^ 2
  } else {
    # param_sdlog = sqrt( -2 * param_meanlog )
    print("Setting lognorm E(X) to 1 | Not ready for calculating SDLOG from MEANLOG")
    res = list(result = -1, err_code = 1, err_msg = "Setting lognorm E(X) to 1 | Not ready for calculating SDLOG from MEANLOG")
    auction__gen_err_msg(res)
  }
  return( list(
    param_meanlog = param_meanlog,
    param_sdlog = param_sdlog
  ) )
}

auction__unobs_dist__std__lognorm <- function(param_meanlog, param_sdlog) {
  # Given meanlog and sdlog, calculate standard deviation of the unobserved heterogeneity
  #   var(x) = ( exp(s ^ 2) - 1 ) * exp(2 * mu + s ^ 2)
  return(sqrt( ( exp(param_sdlog ^ 2) - 1 ) * exp(2 * param_meanlog + param_sdlog ^ 2) ))
}


auction__unobs_dist__exp_val_1__weibull <- function(param_scale = NULL, param_shape = NULL) {
  # Given either scale or shape, set the other such that E(x) = 1
  #   E(x) = scale * GAMMA(1 + 1 / shape)

  if (! is.null(param_shape)) {
    param_scale = 1 / gamma(1 + 1/param_shape)
  } else {
    print("Setting Weibull E(X) to 1 | Not ready for calculating SHAPE from SCALE")
    res = list(result = -1, err_code = 1, err_msg = "Setting Weibull E(X) to 1 | Not ready for calculating SHAPE from SCALE")
    auction__gen_err_msg(res)
  }
  return( list(
    param_scale = param_scale,
    param_shape = param_shape
  ) )
}

auction__unobs_dist__std__weibull <- function(param_scale, param_shape) {
  # Given scale and shape, calculate standard deviation of the unobserved heterogeneity
  #   var(x) = scale ^ 2 * (   GAMMA(1+ 2 / shape)   -   ( gamma(1 + 1 / shape) ) ^ 2    )
  return(sqrt(
    param_scale ^ 2 * ( gamma(1 + 2 / param_shape) - (gamma(1 + 1 / param_shape)) ^ 2 )
  ))
}


auction__unobs_dist__exp_val_1__gamma <- function(param_rate = NULL, param_shape = NULL) {
  # Given either rate or shape, set the other such that E(x) = 1
  #   E(x) = shape / rate

  if (is.null(param_rate)) {
    param_rate = 1 / param_shape
  } else {
    param_shape = 1 / param_rate
  }
  return( list(
    param_rate = param_rate,
    param_shape = param_shape
  ) )
}

auction__unobs_dist__std__gamma <- function(param_rate, param_shape) {
  # Given rate and shape, calculate standard deviation of the unobserved heterogeneity
  #   var(x) = shape / (rate ^ 2)
  return(sqrt( param_shape / (param_rate ^ 2) ))
}



f__ll_parallel__v3 = function(x0, dat_price, dat_num, dat_X, func__prob_distrib, cl){
  # From iteration to iteration, only x0 is changing

  # params = x0
  # v__y = dat_price
  # v__n = dat_num
  # m__h_x = dat_X


  # val__pv_weibull_mu
  # val__pv_weibull_a
  # val__unobs_dist_param
  # val__x_terms

  listIdx = auction__x0_indices()
  # idx__pv_weibull_mu = 1
  # idx__pv_weibull_a = 2
  # idx__unobs_dist_param =3
  # idx__x_terms__start = 4


  # v__mu = x0[idx__pv_weibull_mu]
  # v__alpha = x0[idx__pv_weibull_a]
  # u = x0[idx__unobs_dist_param]


  h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[2]-1)]
  v__h = exp( colSums( h * t(dat_X) ) )

  if (x0[listIdx$unobs_dist_param] <= 0.1)
    return(-Inf) # Check that these hold at estimated values
  else if ( sum (x0[listIdx$pv_weibull_mu] <= 0 ) > 0 )
    return(-Inf)
  else if ( sum( x0[listIdx$pv_weibull_a] <= 0.01 ) > 0)
    return(-Inf)
  else
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/x0[listIdx$pv_weibull_a])
  v__w = dat_price / v__h

  dat = cbind(v__w, dat_num, x0[listIdx$pv_weibull_mu], x0[listIdx$pv_weibull_a], v__gamma_1p1opa)

  # Set E(X) = 1 for UnObserved distribution
  func__prob_distrib = auction__unobs_dist__exp_val_1(
    func__prob_distrib=func__prob_distrib,
    paramVal = x0[listIdx$unobs_dist_param]
    )

  # Run
  v__f_w = parApply(cl = cl,
                    X = dat,
                    MARGIN = 1,
                    FUN = f__funk__v3,
                    func__prob_distrib=func__prob_distrib
                    )

  # Return output
  ### print(paste("v__f_w", -sum(log(v__f_w / v__h))))
  ### print(func__prob_distrib)
  v__f_y = v__f_w / v__h
  return(-sum(log(v__f_y)))
}
f__funk__v3 = function(data_vec, func__prob_distrib){
  val = integrate(vf__w_integrand_z_fast__v3, w_bid=data_vec[1],
                  num_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], func__prob_distrib=func__prob_distrib,
                  lower=0, upper=Inf, abs.tol = 1e-10)
  if(val$message != "OK")
    stop("Integration failed.")
  return(val$value)
}
#vf__w_integrand_z_fast__v3 = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, param_u, func__prob_distrib){
vf__w_integrand_z_fast__v3 = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, func__prob_distrib){
  # Get "x"
  b_z = vf__bid_function_fast__v3(price=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z

  # Add "x"
  ### pdf_list_entry$argList$x = u_z
  if ( length(func__prob_distrib$argList)==1 && is.null(func__prob_distrib$argList[[1]]) ) {
    func__prob_distrib$argList = list(x = u_z)
  } else {
    func__prob_distrib$argList$x = u_z
  }

  #Run
  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    do.call(
      match.fun(func__prob_distrib$funcName),
      func__prob_distrib$argList
    )
  ### dlnorm(u_z, meanlog=(-param_u^2*1/2), sdlog = param_u) # Note: can swap for different distributions

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}

vf__bid_function_fast__v3 = Vectorize(FUN = f__bid_function_fast__v3,vectorize.args = "price")

f__bid_function_fast__v3 = function(price, num_bids, mu, alpha, gamma_1p1oa){

  if (exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*price)^alpha) == 0) {
    return(price + mu/alpha*(num_bids-1)^(-1/alpha)*1/gamma_1p1oa*
             ((num_bids-1)*(gamma_1p1oa/mu*price)^alpha)^(1/alpha-1))
  }

  price + 1/alpha*(mu/gamma_1p1oa)*(num_bids-1)^(-1/alpha)*
    pgamma((num_bids-1)*(1/(mu/gamma_1p1oa)*price)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*price)^alpha)
  # Check gamma(1/alpha) part
}



auction__x0_indices <- function() {
  return( list(
    pv_weibull_mu = 1,
    pv_weibull_a = 2,
    unobs_dist_param = 3,
    x_terms__start = 4
    ) )
}
auction__x0_stepsizes <- function() {
  return( list(
    pv_weibull_mu = 1,
    pv_weibull_a = 0.1,
    unobs_dist_param = 1,
    x_terms__start = 0.1,
    x_terms__other = 1
    ) )
}



listInputPDF = list(dlnorm=list()) # working
listInputPDF = list(dweibull=list()) # working
listInputPDF = NULL

auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 2,
           func_list__unobs_distrib = listInputPDF
           )



































auction <- function (input_data = NULL, initial_guess = NULL,
                           generate_data = FALSE, pdf_list = NULL,
                           num_cores = 1) {
  # INPUTS
  #   input_data
  #
  #
  #   initial_guess [list]
  #
  #
  #   generate_data [boolean]
  #
  #
  #   pdf_list [list]
  #     define the probability density function(s) to use
  #     input must be list of the form
  #       pdf_list = list( pdf1 = list( param_a = , param_b = ), pdf2 = list( ... ) )
  #
  #     currently supported PDFs are:
  #       dnorm - normal distribution
  #       dlnorm - log normal distribution
  #
  #     valid parameters for each PDF may be found at
  #       https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Distributions.html
  #
  #
  #   num_cores [int]
  #     Number of parallel workers to spawn.
  #
  #
  #
  # OUTPUTS
  #   results
  #     Your results
  #
  #   err_code
  #     0 = no error
  #     1 = unable to load one or more required packages
  #     2 = one or more input parameters is invalid
  #     3 = Input data or initial guess were not given
  #     4 = Input data is not correctly structured
  #
  #   err_msg
  #     Message accompanying the error code


  # Validate environment and inputs
  #   Initialize
  res <- list(result = -1, err_code = 0, err_msg = "")

  #   Validate environment
  #     Load required libraries
  res = auction__load_packages(res)

  # Check inputs
  #   Check parameter data types for non-null default value
  res <- auction__check_input__generate_data(generate_data, res)
  res <- auction__check_input__num_cores(num_cores, res)
  #   Check PDFs specified and untangle pdf_list and res
  #     (1) run check
  res <- auction__check_input__pdf_list(pdf_list, res)
  #     (2) untangle
  pdf_list <- res$inp
  res <- res$res


  # Check to see if we should error out
  if ( res['err_code'] != 0 ) {
    return(res)
  }

  #   Check input data and initial guess
  if ( ! generate_data ) {
    if ( is.null(input_data) ) {
      res['err_code'] = 3
      res['err_msg'] = "No data given. Use 'generate_data=TRUE' to generate sample input data"
    }
    if ( is.null(initial_guess) ) {
      res['err_code'] = 3
      res['err_msg'] = "No initial guess given"
    }
    if ( res['err_code'] != 0 ) {
      return(res)
    }
  }



  # Prepare data
  #   Part 1) Inspect data
  if (generate_data) {
    dat = generate__data()
  } else {

    dat = input_data

    res = auction__check_input__input_data__columns(dat, res)
    if ( res['err_code'] != 0 ) {
      return(res)
    }
  }

  #   Part 2) Gather components
  v__y = dat$v__y
  v__n = dat$v__n
  m__h_x = dat$m__h_x


  # Get initial guess for convergence
  if ( generate_data ) {
    x0 = genereate__initial_guess(v__y, v__n, m__h_x)
  } else {
    x0 = initial_guess
  }




  # Set up parallelization
  cl = makeCluster(num_cores)
  clusterExport(cl,varlist=c("vf__bid_function_fast",
                             "vf__w_integrand_z_fast__2",
                             "f__funk"))


  #   Get convergence conditions
  optim_control = f_conv_conditions(x0)


  print(optim_control)
  return(NULL)


  # Run
  # run_result = optim(par=x0, fn=f__ll_parallel, control=optim_control,
  #                y=v__y, n=v__n, h_x=m__h_x, cl=cl)

  run_result = list()
  for (list_entry in names(pdf_list)) {
    # Initialize results
    run_result[[as.character(list_entry)]] = NULL

    # Prepare inputs
    pdf_list_entry = list(funcName=as.character(list_entry), argList=pdf_list[[as.character(list_entry)]])

    # Run
    print(paste("Running |", as.character(list_entry)))

    run_result[[as.character(list_entry)]] = optim(par=x0, fn=f__ll_parallel__2, control=optim_control,
                       y=v__y, n=v__n, h_x=m__h_x,
                       pdf_list_entry=pdf_list_entry,
                       cl=cl)

    print(paste("End of run |", as.character(list_entry)))
  }


  # Release resources
  stopCluster(cl)


  # Inspect result

  # If inspection is fine, add to res
  res$result = run_result


  # Return result
  return(res)
}

#' @param obs, number of observations to draw,
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
generate__data <- function(obs = 200) {
  # For testing purposes, we will generate sample data

  # ensure that obs is integer greater than 0
  set.seed(301)
  # data = # Generate some data
  # y, n, x1, x2: positive
  # n: discrete and > 1
  # y is some function of n, x1, x2

  #obs = 200
  obs = 20
  w = rlnorm(obs)
  x1 = rlnorm(obs) + 0.5*w
  x2 = 0.1*rlnorm(obs) + 0.3*w
  e = 2*rlnorm(obs)
  n = sample(2:10, obs, replace=TRUE)
  y = 10 - 0.5*n + x1 + x2 + e
  data = data.frame(cbind(y, n, x1, x2))
  #plot(n, y)

  v__y = data$y
  v__n = data$n
  m__h_x = as.matrix(cbind(log(data$x1),log(data$x2)))

  return(
    list(
      v__y = v__y, v__n = v__n, m__h_x = m__h_x
      )
  )
}


###########################################################################
#
###########################################################################



auction__load_packages <- function (res) {
  if ( res['err_code'] == 0 ) {
    listRequiredPackages = c('parallel')
    listMissingPackages = c()
    for (reqPkg in listRequiredPackages) {
      if ( ! require(reqPkg, character.only=TRUE)) {
        listMissingPackages = c(listMissingPackages, reqPkg)
      }
    }
    if ( length(listMissingPackages) > 0 ) {
      res['err_code'] = 1
      res['err_msg'] = paste0("Unable to load the following packages: ",
                              paste(listMissingPackages, collapse=','))
    }
  }
  return(res)
}

auction__check_input__generate_data <- function(inp, res) {
  if ( res['err_code'] == 0 ) {
    if ( ! is.logical(inp) ) {
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'generate_data' | Must be boolean"
    }
  }
  return(res)
}
auction__check_input__num_cores <- function(inp, res) {
  if ( res['err_code'] == 0 ) {
    if ( ( ! is.numeric(inp) ) || ( inp%%1 != 0) ||  ( inp < 0 ) ) {
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'num_cores' | Must be natural number"
    } else {
      # Check cores available
      num_cores_detected = detectCores(res)
      if ( inp > num_cores_detected ) {
        print(paste0("Warning: You have requested ", inp,
                     " cores but only have ", num_cores_detected,
                     " cores available"))
      }
    }
  }
  return(res)
}
auction__check_input__input_data__columns <- function(inp, res) {
  if ( res['err_code'] == 0 ) {
    colList = names(inp)
    reqList = c("v__y", "v__n", "m__h_x")
    missingList = c()
    for (reqCol in reqList) {
      if ( ! reqCol %in% colList ) {
        missingList = c(missingList, reqCol)
      }
    }
    if ( length(missingList) > 0 ) {
      res['err_code'] = 4
      res['err_msg'] = paste0("Missing the following columns: ",
                              paste(missingList, collapse=','))
    }
  }
  return(res)
}

auction__check_input__pdf_list <- function(inp, res) {
  if ( res['err_code'] == 0 ) {
    if (is.null(inp)) {
      inp = list(dlnorm=NULL)
    }

    if (! is.list(inp) || length(inp) == 0) {
      res['err_code'] = 2
      res['err_msg'] = "Invalid input for 'pdf_list' | Must be NULL for default PDF or list as per documentation"
    } else {
      pdf_list = inp

      # Get list of valid PDFs
      valid_pdf_list = auction__valid_opt__pdf_list()

      # Remove PDFs that are not valid
      pdf_list = pdf_list[names(pdf_list) %in% names(valid_pdf_list)]

      if (length(pdf_list)==0) {
        res['err_code'] = 2
        res['err_msg'] = "Invalid PDFs specified for 'pdf_list'"
      } else {

        # Check to see if parameter requirements are met
        for (funcName in names(pdf_list)) {
          I_validPDF = TRUE

          curPDF = pdf_list[[as.character(funcName)]]

          # Remove parameters that are not valid
          if (is.list(curPDF)) {
            # Remove by parameter name
            curPDF = curPDF[names(curPDF) %in% names(valid_pdf_list[[as.character(funcName)]])]

            if (length(names(curPDF)) != 0){
              # Remove by parameter datatype
              for (paramName in names(curPDF)) {
                # Get parameter value being tested
                paramVal = curPDF[[as.character(paramName)]]
                #   If NULL, then remove
                #   Otherwise, if not boolean nor numeric, then wrap in quotes
                if (is.null(paramVal)) {
                  # Remove and continue
                  curPDF = curPDF[names(curPDF) != paramName]
                  next
                } else if ((! is.logical(paramVal)) && (! is.numeric(paramVal))) {
                  # Wrap in quotes
                  paramVal = paste("'", paramVal, "'", sep="")
                }
                # Get required parameter data type
                paramReqType = valid_pdf_list[[as.character(funcName)]][[as.character(paramName)]][['type']]

                # Check
                checkStr = paste("is.", paramReqType, "(", paramVal, ")", sep="")
                I_check = eval(parse(text=checkStr))

                if (! I_check) {
                  # Remove
                  curPDF = curPDF[names(curPDF) != paramName]
                }

              }
            }
          }

          # Check for required parameters
          #   curPDF or length(curPDF)==0 is okay if no required parameters
          if (is.null(curPDF) || (is.list(curPDF))) {
            # Look for requirements
            req_list = c()
            validPDF = valid_pdf_list[[as.character(funcName)]]
            if (is.list(validPDF) && length(validPDF) > 0) {
              for (paramName in names(validPDF)) {

                if (validPDF[[as.character(paramName)]][['req']]) {
                  req_list = c(req_list, paramName)
                }

              }
            }

            # Check to see if requirements are met
            if (length(req_list) > 0) {
              if (is.null(curPDF) || (length(names(curPDF))==0)) {
                I_validPDF = FALSE
              } else {

                I_req = reqList %in% names(curPDF)
                if (! all(I_req)) {
                  I_validPDF = FALSE
                  print(paste(
                    "Ignoring PDF ", funcName, " because missing required parameters: ",
                    paste(reqList[! I_req], collapse=','), sep = ""
                    ))
                }

              }
            }

          } else {
            I_validPDF = FALSE
          }

          if (! I_validPDF) {
            pdf_list = pdf_list[names(pdf_list) != funcName]
          } else {
            if (length(names(curPDF))==0) {
              curPDF = list(NULL)
            }

            pdf_list[[as.character(funcName)]] = curPDF

          }
        }

        if (length(pdf_list)==0) {
          res['err_code'] = 2
          res['err_msg'] = "No PDFs specified in 'pdf_list' have valid required parameters"
        }

      }

      inp = pdf_list

    }
  }
  return(list(res=res,inp=inp))
}
auction__valid_opt__pdf_list <- function() {
  # Parameters have two values, (1) whether they are required and (2) what datatype must they be
  #   If a parameter is not required, than (1) is NULL
  return(list(
    dcauchy=list(
      location=list(type='numeric',req=FALSE),
      scale=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      ),
    dchisq=list(
      df=list(type='numeric',req=TRUE),
      ncp=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      ),
    dexp=list(
      rate=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      ),
    dgamma=list(
      shape=list(type='numeric',req=TRUE),
      rate=list(type='numeric',req=FALSE),
      scale=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      ),
    dlnorm=list(
      meanlog=list(type='numeric',req=FALSE),
      sdlog=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      ),
    dnorm=list(
      mean=list(type='numeric',req=FALSE),
      sd=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      ),
    dweibull=list(
      shape=list(type='numeric',req=TRUE),
      scale=list(type='numeric',req=FALSE),
      log=list(type='logical',req=FALSE)
      )
    ))
}


# auction__check_input <- function(inp, inpName, inpTypeList, res) {
#   if ( res['err_code'] == 0 ) {
#     errNum = 2
#     I_valid = FALSE
#     for ( inpType in inpTypeList ) {
#       inpTypeSplit = strsplit(inpType, "_")[[1]]
#
#       if ( inpTypeSplit[1] == 'numeric' ) {
#
#         if ( length(inpTypeSplit) == 1) {
#           if ( ! is.numeric(inp) ) {
#             res['err_code'] = errNum
#           }
#         } else {
#           if ( inpTypeSplit[2] == 'natural' ) {
#             if ( ( ! is.numeric(inp) ) || ( inp%%1 != 0) ||  ( inp < 0 ) ) {
#               res['err_code'] = errNum
#             }
#           }
#         }
#
#       } else if ( inpTypeSplit[1] == 'string' ) {
#
#       }
#
#     }
#
#
#     # if ( ( ! is.numeric(inp) ) || ( inp%%1 != 0) ||  ( inp < 0 ) ) {
#     #   return(list(
#     #     result=-1,
#     #     err_code=2,
#     #     err_msg="Invalid input for 'num_cores'"
#     #   ))
#     # }
#     if ( res['err_code'] == 2 ) {
#       res['err_msg'] = paste0("Invalid input for '", inpName, "'")
#     }
#   }
#   return(res)
# }






genereate__initial_guess <- function(v__y, v__n, m__h_x) {
  # Find starting point / intial guess
  #
  # Input:
  #   v__y
  #
  #   v__n
  #
  #   m__h_x
  #
  # Output:
  #   Returned list is of the form
  #     (mu, a, u, h[1], h[2], ... )

  # Currently, fix initial guess
  mu_init = 8
  a_init = 2
  u_init = 0.5
  h_init = rep(0.5, dim(m__h_x)[2])
  h_init[1] = 0.4
  h_init[2] = 0.5

  return(c(mu_init, a_init, u_init, h_init))
}



f_conv_conditions <- function(x_initialGuess) {
  # Define convergence conditions
  #
  # Input:
  #   x_initialGuess
  #     must be a numeric array of at least length 5
  return(
    list(
      maxit = 2000,
      parscale = c(1, 0.1, 1, 0.1, rep(1, length(x_initialGuess) - 4) ))
    )
}



###########################################################################
# Estimation Functions
###########################################################################
f__bid_function_fast = function(cost, num_bids, mu, alpha, gamma_1p1oa){

  if (exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha) == 0) {
    return(cost + mu/alpha*(num_bids-1)^(-1/alpha)*1/gamma_1p1oa*
             ((num_bids-1)*(gamma_1p1oa/mu*cost)^alpha)^(1/alpha-1))
  }

  cost + 1/alpha*(mu/gamma_1p1oa)*(num_bids-1)^(-1/alpha)*
    pgamma((num_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha, 1/alpha, lower=FALSE)*
    gamma(1/alpha)*
    1/exp(-(num_bids-1)*(1/(mu/gamma_1p1oa)*cost)^alpha)
  # Check gamma(1/alpha) part
}

vf__bid_function_fast = Vectorize(FUN = f__bid_function_fast,vectorize.args = "cost")

vf__w_integrand_z_fast = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, param_u){

  b_z = vf__bid_function_fast(cost=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z

  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    dlnorm(u_z, meanlog=(-param_u^2*1/2), sdlog = param_u) # Note: can swap for different distributions

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}

vf__w_integrand_z_fast__2 = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, param_u, pdf_list_entry){
  # Get "x"
  b_z = vf__bid_function_fast(cost=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z

  # Add "x"
  ### pdf_list_entry$argList$x = u_z
  if ( length(pdf_list_entry$argList)==1 && is.null(pdf_list_entry$argList[[1]]) ) {
    pdf_list_entry$argList = list(x = u_z)
  } else {
    pdf_list_entry$argList$x = u_z
  }

  #Run
  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    do.call(
      match.fun(pdf_list_entry$funcName),
      pdf_list_entry$argList
    )
    ### dlnorm(u_z, meanlog=(-param_u^2*1/2), sdlog = param_u) # Note: can swap for different distributions

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)
}



f__funk = function(data_vec, param_u, pdf_list_entry){

  # val = integrate(vf__w_integrand_z_fast, w_bid=data_vec[1],
  val = integrate(vf__w_integrand_z_fast__2, w_bid=data_vec[1],
                  num_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], param_u=param_u, pdf_list_entry=pdf_list_entry,
                  lower=0, upper=Inf, abs.tol = 1e-10)
  #                   rel.tol = .Machine$double.eps^0.7)
  if(val$message != "OK")
    stop("Integration failed.")
  return(val$value)
}



f__ll_parallel = function(x0, y, n, h_x, cl){
  params = x0
  v__y = y
  v__n = n
  m__h_x = h_x

  v__mu = params[1]
  v__alpha = params[2]
  u = params[3]

  h = params[4:( 3 + dim(m__h_x)[2] )]
  v__h = exp( colSums( h * t(m__h_x) ) )

  if (u <= 0.1)
    return(-Inf) # Check that these hold at estimated values
  else if ( sum (v__mu <= 0 ) > 0 )
    return(-Inf)
  else if ( sum( v__alpha <= 0.01 ) > 0)
    return(-Inf)
  else
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/v__alpha)
    v__w = v__y / v__h
    dat = cbind(v__w, v__n, v__mu, v__alpha, v__gamma_1p1opa)

    v__f_w = parApply(cl = cl, X = dat, MARGIN = 1, FUN = f__funk, param_u = u)
    v__f_y = v__f_w / v__h
    return(-sum(log(v__f_y)))
}
f__ll_parallel__2 = function(x0, y, n, h_x, pdf_list_entry, cl){
  # From iteration to iteration, only x0 is changing
  params = x0
  v__y = y
  v__n = n
  m__h_x = h_x

  v__mu = params[1]
  v__alpha = params[2]
  u = params[3]

  h = params[4:( 3 + dim(m__h_x)[2] )]
  v__h = exp( colSums( h * t(m__h_x) ) )

  if (u <= 0.1)
    return(-Inf) # Check that these hold at estimated values
  else if ( sum (v__mu <= 0 ) > 0 )
    return(-Inf)
  else if ( sum( v__alpha <= 0.01 ) > 0)
    return(-Inf)
  else
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/v__alpha)
  v__w = v__y / v__h


  dat = cbind(v__w, v__n, v__mu, v__alpha, v__gamma_1p1opa)

  # print("dat")
  # print(dat)

  v__f_w = parApply(cl = cl, X = dat, MARGIN = 1,
                    FUN = f__funk,
                    param_u = u,
                    pdf_list_entry=pdf_list_entry
  )

  print("v__f_w")
  # print(v__f_w)
  print(-sum(log(v__f_w / v__h)))

  v__f_y = v__f_w / v__h
  return(-sum(log(v__f_y)))
}







###########################################################################
# Run code
###########################################################################

# listInputPDF = list(
#   fakePDF=list(fake1=1, fake2=2),
#   dlnorm=list(sdlog=1)
#   ) # working
# listInputPDF = NULL # working -> defaults to dlnorm
# listInputPDF = list(dlnorm=NULL) # working
listInputPDF = list(dlnorm=list()) # working
# listInputPDF = list(dlnorm=NULL, dnorm=NULL) # working
# listInputPDF = list(dgamma=NULL) # working

res = auction(generate_data=TRUE, num_cores=2, pdf_list=listInputPDF)
print(res)
