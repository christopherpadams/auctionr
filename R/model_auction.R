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
  res = auction__init_env(num_cores=num_cores)
  #   Unpack 'res'
  num_cores = res$num_cores
  res = res$res

  # Validate that data, winning_bid and number_of_bids all match
  res = auction__input_data_consistency(res=res,
                                        dat=dat,
                                        colName_price=winning_bid,
                                        colName_num=number_of_bids)

  # Check consistency of initial guess and distribution information for unobserved heterogeneity
  res = auction__check_input__unobs_distrib_AND_initial_guess(func_list=func_list__unobs_distrib,
                                                              initial_guess=initial_guess,
                                                              dat=dat,
                                                              res=res)
  #     unpack results
  func_list__unobs_distrib = res$func_list
  initial_guess = res$initial_guess
  res = res$res

  # Set up parallelization
  cl = makeCluster(num_cores)
  clusterExport(cl,varlist=c("vf__bid_function_fast__v3",
                             "vf__w_integrand_z_fast__v3",
                             "f__funk__v3"))
  # Run
  run_result = list()
  for (funcName in names(func_list__unobs_distrib)) {
    sFuncName = as.character(funcName)

    # Initialize results
    run_result[[sFuncName]] = NULL

    # Prepare inputs
    #   unobserved heterogeneity distribution
    func__unobs_distrib = func_list__unobs_distrib[[sFuncName]]

    #   initial guess, step sizes, max iterations
    conv_params = auction__conv__get_params(
      func__unobs_distrib=func_list__unobs_distrib[[sFuncName]],
      initial_guess=initial_guess )
    #     unpack
    I_unobs__const = conv_params$I_unobs__const
    x0 = conv_params$initial_guess
    conv_ctrl = conv_params$conv_ctrl

    #   set I_unobs__const
    func__unobs_distrib$I_unobs__const = I_unobs__const
    func_list__unobs_distrib[[sFuncName]]$I_unobs__const = I_unobs__const

    # Run
    print(paste("Running |", sFuncName))

    run_result[[sFuncName]] = optim(par=x0,
                                    fn=f__ll_parallel__v3,
                                    control=conv_ctrl,
                                    dat_price=dat[[winning_bid]],
                                    dat_num=dat[[number_of_bids]],
                                    dat_X=dat[['x_terms']],
                                    func__prob_distrib=func__unobs_distrib,
                                    cl=cl)

    print(paste("End of run |", sFuncName))
    print(run_result[[sFuncName]])
  }
  # Release resources
  stopCluster(cl)
  # Inspect result
  res$result = auction__output_org(run_result=run_result, func_list__unobs_distrib=func_list__unobs_distrib)
  # Return result
  return(res)
}


# Called by auction_v2
auction__output_org <- function(run_result, func_list__unobs_distrib) {

  # Get standard deviations of the unobserved heterogeneity
  run_result__raw = list()
  run_result__df = list()
  for (funcName in names(run_result)) {
    sFuncName = as.character(funcName)

    # What was the parameter being changed?
    paramName = func_list__unobs_distrib[[sFuncName]]$argName_ctrl
    # What was its final value?
    if (func_list__unobs_distrib[[sFuncName]]$I_unobs__const) {
      paramVal = func_list__unobs_distrib[[sFuncName]]$argList[[paramName]]
    } else {
      listIdx = auction__x0_indices()
      paramVal = run_result[[sFuncName]]$par[listIdx$unobs_dist_param]
    }

    # Use auction__unobs_dist__exp_val_1() to get other parameter plus its value
    func__prob_distrib = auction__unobs_dist__exp_val_1(func__prob_distrib = func_list__unobs_distrib[[sFuncName]],
                                                        paramVal = paramVal)

    # Get standard deviation of unobserved heterogeneity
    unobs_std = auction__unobs_dist__std(func__prob_distrib = func__prob_distrib)

    # Build output list
    run_result__raw[[sFuncName]] = list()
    run_result__raw[[sFuncName]]$invLogLikelyhood = run_result[[sFuncName]]$value
    run_result__raw[[sFuncName]]$funcName = func__prob_distrib$funcName
    run_result__raw[[sFuncName]]$argList = func__prob_distrib$argList
    run_result__raw[[sFuncName]]$argName_ctrl = func__prob_distrib$argName_ctrl
    run_result__raw[[sFuncName]]$I_unobs__const = func__prob_distrib$I_unobs__const
    run_result__raw[[sFuncName]]$std = unobs_std

    # Build dataframe result to organize by StdDev
    run_result__df[[sFuncName]] = run_result__raw[[sFuncName]]
    #   Remove parameters which are lists themselves
    run_result__df[[sFuncName]]$argList = NULL
    ##  (is.list() ) not working, try length()>1 ?
    # for (argKey in names(run_result__df[[sFuncName]])) {
    #   if (is.list( run_result__df[[sFuncName]][[as.character(argKey)]] )) {
    #     run_result__df[[sFuncName]][[as.character(argKey)]] = NULL
    #   }
    # }

  }
  # Convert list to dataframe
  run_result__df = do.call(rbind.data.frame, run_result__df)
  # Sort by standard deviation
  run_result__df = run_result__df[with(run_result__df, order(std)), ]
  # Add distribution parameters back in
  run_result__df['funcParam1'] = NaN
  run_result__df['funcParam1_val'] = NaN
  run_result__df['funcParam2'] = NaN
  run_result__df['funcParam2_val'] = NaN
  for (funcName in names(run_result)) {
    sFuncName = as.character(funcName)

    run_result__df[sFuncName, 'funcParam1'] = run_result__raw[[sFuncName]]$argName_ctrl

    run_result__df[sFuncName, 'funcParam1_val'] = run_result__raw[[sFuncName]]$argList[[
      run_result__raw[[sFuncName]]$argName_ctrl ]]

    run_result__df[sFuncName, 'funcParam2'] = names(run_result__raw[[sFuncName]]$argList)[
      ! names(run_result__raw[[sFuncName]]$argList) %in% run_result__df[sFuncName, 'funcParam1'] ]

    run_result__df[sFuncName, 'funcParam2_val'] =
      run_result__raw[[sFuncName]]$argList[[run_result__df[sFuncName, 'funcParam2']]]
  }
  # Re-org columns
  #   [funcName, std, invLogLikelyhood, etc]
  listColNames = c('funcName', 'std', 'invLogLikelyhood', 'argName_ctrl')
  run_result__df = run_result__df[c(listColNames, colnames(run_result__df)[! colnames(run_result__df) %in% listColNames])]

  return(run_result__df)
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

auction__init_env <- function (num_cores) {
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

  return( list(
    res = res,
    num_cores = num_cores
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
    num_cores__avail = detectCores()
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

auction__check_input__unobs_distrib_AND_initial_guess <- function(func_list, initial_guess, dat, res) {

  # Prepare 'func_list__unobs_distrib'
  #   Initial prep
  res = auction__check_input__unobs_distrib__prep(func_list=func_list,
                                                        res=res)
  #       unpack
  func_list = res$func_list
  res = res$res
  #   Clean-up 'func_list__unobs_distrib'
  #     Remove invalid distributions
  #     Remove invalid parameters
  res = auction__check_input__unobs_distrib__validate_func(func_list=func_list, res=res)
  #       unpack
  func_list = res$func_list
  res = res$res
  #   Make sure a control parameter is set
  res = auction__check_input__unobs_distrib__set_ctrl_arg(func_list=func_list, res=res)
  #       unpack
  func_list = res$func_list
  res = res$res

  # Prepare 'initial_guess'
  #   separate into initial guess values for parameters NOT related to
  #   unobserved heterogeneity distributions, and the parameters which
  #   are related
  nParams_dat=auction__get_num_columns__dat(dat)
  initial_guess = auction__check_input__inital_guess__prep(initial_guess=initial_guess,
                                                           nParams_dat=nParams_dat,
                                                           res=res)
  #   unpack
  initial_guess__base = initial_guess$base
  initial_guess__unobs =initial_guess$unobs

  # Make sure we have an initial guess for all unobserved heterogeneity distributions
  func_list = auction__unobs_distrib__set_init_guess(func_list=func_list,
                                                                  initial_guess=initial_guess__unobs)

  # Set step sizes for unobserved heterogeneity distributions
  func_list = auction__unobs_distrib__set_step_size(func_list)

  # Clean-up
  #   Reset 'initial_guess' to 'initial_guess__base'
  initial_guess = initial_guess__base

  return( list(
    res = res,
    func_list__unobs_distrib = func_list,
    initial_guess = initial_guess
    ) )
}

auction__check_input__unobs_distrib__prep <- function(func_list, res) {

  if ( !is.null(func_list) &&
       !is.vector(func_list) &&
       !is.list(func_list) ) {
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'func_list__unobs_distrib' | Must be NULL, a string vector, or a list"
    auction__gen_err_msg(res)
  }

  if (! is.list(func_list) && ! is.character(func_list)) {
    # No distributions given, set to default
    func_list = list( Default=auction__check_input__unobs_distrib__default_func() )
  } else if (is.character(func_list)) {
    # Distributions given as strings
    tmp = list()
    for (sFuncName in func_list) {
      tmp[[sFuncName]] = auction__check_input__unobs_distrib__default_func(sFuncName)
    }
    func_list = tmp
  } else if (is.list(func_list)) {
    # Received a list, make sure it has appropriate fields
    if (length(func_list)==0) {
      # No distributions given, set to default
      #   log normal
      func_list = list( Default=auction__check_input__unobs_distrib__default_func() )
    } else {
      for (funcName in names(func_list)) {
        sFuncName = as.character(funcName)

        if (! is.list(func_list[[sFuncName]]) || length(func_list[[sFuncName]])==0) {
          # Not sure what this would be instead
          print(paste("'func_list__unobs_distrib' key '", sFuncName, "' has unknown value. Building default entry", sep=''))
          func_list[[sFuncName]] = auction__check_input__unobs_distrib__default_func(sFuncName)
        } else {
          # Make sure we have all parameters

          if (! 'funcName' %in% names(func_list[[sFuncName]])) {
            func_list[[sFuncName]]$funcName = sFuncName
          }

          func_list[[sFuncName]] = auction__check_input__unobs_distrib__default_func(func_list[[sFuncName]])
        }
      }
    }

  }

  # # Make sure a control parameter specified
  # for (funcName in names(func_list)) {
  #   sFuncName = as.character(funcName)
  #
  #   paramName_ctrl = NULL
  #   if (! is.nan(func_list[[sFuncName]]$argName_ctrl) &&
  #       is.character(func_list[[sFuncName]]$argName_ctrl) &&
  #       length(func_list[[sFuncName]]$argName_ctrl)==1) {
  #
  #     paramName_ctrl = func_list[[sFuncName]]$argName_ctrl
  #
  #     # Does the value exist within $argList?
  #     if (is.list(func_list[[sFuncName]]$argList) && length(func_list[[sFuncName]]$argList) > 0) {
  #
  #     }
  #
  #   } else {
  #     # Do we have any entries within $arglist?
  #
  #   }
  # }

  return(list(func_list=func_list, res=res))
}
auction__check_input__unobs_distrib__default_func <- function(func__unobs_distrib = NULL) {
  if (is.null(func__unobs_distrib) || length(func__unobs_distrib) == 0) {
    func__unobs_distrib = list( funcName='dlnorm',
                                argList=list(),
                                argName_ctrl=NaN,
                                initial_guess=NaN,
                                step_size=NaN )
  } else if (is.character(func__unobs_distrib) && length(func__unobs_distrib)==1) {
    func__unobs_distrib = list( funcName=func__unobs_distrib,
                                argList=list(),
                                argName_ctrl=NaN,
                                initial_guess=NaN,
                                step_size=NaN )
  } else if (is.list(func__unobs_distrib) && length(func__unobs_distrib) > 0) {
    # Make sure we have all parameters
    listKeys = names(func__unobs_distrib)

    if (! 'funcName' %in% listKeys) {
      print("auction__check_input__unobs_distrib__default_func | $funcName not defined! Resetting to default")
      func__unobs_distrib = list( funcName='dlnorm',
                                  argList=list(),
                                  argName_ctrl=NaN,
                                  initial_guess=NaN,
                                  step_size=NaN )
    } else {
      if (! 'argList' %in% listKeys) {
        func__unobs_distrib$argList = list()
      } else if (is.character(func__unobs_distrib$argList)) {
        # Convert to list
        tmp = func__unobs_distrib$argList
        func__unobs_distrib$argList = list()
        for (argKey in tmp) {
          func__unobs_distrib$argList[[argKey]] = NaN
        }
      }
      if (! 'argName_ctrl' %in% listKeys) {
        func__unobs_distrib$argName_ctrl = NaN
      }
      if (! 'initial_guess' %in% listKeys) {
        func__unobs_distrib$initial_guess = NaN
      }
      if (! 'step_size' %in% listKeys) {
        func__unobs_distrib$step_size = NaN
      }
    }

  } else {
    print("auction__check_input__unobs_distrib__default_func | Unknown entry for 'func__unobs_distrib' | Setting to default")
    func__unobs_distrib = list( funcName='dlnorm',
                                argList=list(),
                                argName_ctrl=NaN,
                                initial_guess=NaN,
                                step_size=NaN )
  }
  return(func__unobs_distrib)
}
auction__check_input__unobs_distrib__validate_func <- function(func_list, res) {
  # Check distributions given in 'func_list'
  #   Remove based on function name
  #     Get list of valid function
  func_list__valid = auction__valid_opt__func_list()
  listValid = names(func_list__valid)
  #     Remove non-whitelisted functions
  listRemove = c()
  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)
    if (! any(listValid %in% func_list[[sFuncName]]$funcName)) {
      listRemove = c(listRemove, sFuncName)
    }
  }
  if (length(listRemove) > 0) {
    for (funcName in listRemove) {
      func_list = func_list[! names(func_list) %in% funcName]
    }
  }

  #   For remaining distributions, remove any parameters that are not whitelisted
  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)

    if (is.list(func_list[[sFuncName]]$argList) && length(names(func_list[[sFuncName]]$argList))>0) {
      func_list[[sFuncName]]$argList = func_list[[sFuncName]]$argList[
        names(func_list[[sFuncName]]$argList) %in% names(func_list__valid[[func_list[[sFuncName]]$funcName]])
      ]

      if (length(func_list[[sFuncName]]$argList)==0) {
        func_list[[sFuncName]]$argList = list()
      }
    }
  }

  #   For remaining distributions, remove any whose parameters are the wrong datatype
  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)

    if (is.list(func_list[[sFuncName]]$argList) && length(names(func_list[[sFuncName]]$argList))>0) {

      listRemove = c()
      for (argKey in names(func_list[[sFuncName]]$argList)) {
        # Check if not NaN
        if (! is.nan(func_list[[sFuncName]]$argList[[argKey]])) {
          # Get required parameter data type
          argType_req = func_list__valid[[func_list[[sFuncName]]$funcName]][[as.character(argKey)]][['type']]
          # Check
          checkStr = paste("is.", argType_req, "(", func_list[[sFuncName]]$argList[[argKey]], ")", sep="")
          I_check = eval(parse(text=checkStr))
          if (! I_check) {
            listRemove = c(listRemove, argKey)
          }
        }
      }

      if (length(listRemove) > 0) {
        for (argKey in listRemove) {
          func_list[[sFuncName]]$argList =
            func_list[[sFuncName]]$argList[ ! names(func_list[[sFuncName]]$argList) %in% argKey ]
        }

        if (length(func_list[[sFuncName]]$argList)==0) {
          func_list[[sFuncName]]$argList = list()
        }
      }
    }
  }

  #   For remaining distributions, remove any functions that are missing required parameters
  listRemove = c()
  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)

    # Check if this function has any required parameters
    listReq = c()
    for (paramName in names(func_list__valid[[func_list[[sFuncName]]$funcName]])) {
      if (func_list__valid[[func_list[[sFuncName]]$funcName]][[paramName]][['req']]) {
        listReq = c(listReq, paramName)
      }
    }
    if (length(listReq)>0) {
      for (paramName_req in listReq) {
        if (! names(func_list[[sFuncName]]$argList) %in% paramName_req ) {
          listRemove = c(listRemove, sFuncName)
          print(paste("Removing '", sFuncName, "' because you did not provide required parameter(s)", sep=''))
          break
        } else if (is.nan(func_list[[sFuncName]]$argList[[paramName_req]]) ||
                   length(func_list[[sFuncName]]$argList[[paramName_req]]==0)) {
          listRemove = c(listRemove, sFuncName)
          print(paste("Removing '", sFuncName, "' because no valid value found for required parameter(s)", sep=''))
          break
        }
      }
    }
  }
  if (length(listRemove) > 0) {
    for (funcName in listRemove) {
      func_list = func_list[! names(func_list) %in% funcName]
    }
  }

  return( list( func_list=func_list, res=res ) )
}
auction__check_input__unobs_distrib__set_ctrl_arg <- function(func_list, res) {

  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)

    # Initialize
    paramName_ctrl = NULL
    paramName_ctrl__default = auction__unobs_dist__ctrl_param__default(func_list[[sFuncName]]$funcName)

    # Search for existing entry
    if (! is.nan(func_list[[sFuncName]]$argName_ctrl) &&
        is.character(func_list[[sFuncName]]$argName_ctrl) &&
        length(func_list[[sFuncName]]$argName_ctrl)==1) {
      paramName_ctrl = func_list[[sFuncName]]$argName_ctrl
    }
    if (is.null(paramName_ctrl)) {
      # Search within $argList
      if (is.list(func_list[[sFuncName]]$argList) && length(func_list[[sFuncName]]$argList) > 0) {
        if (length(func_list[[sFuncName]]$argList)==1) {
          # Only one param listed, it must be the control param
          paramName_ctrl = names(func_list[[sFuncName]]$argList)
        } else {
          # More than one param listed, do any have values?
          listKeys_arg = names(func_list[[sFuncName]]$argList)

          if (paramName_ctrl__default %in% listKeys_arg &&
              length(func_list[[sFuncName]]$argList[[paramName_pref]]) != 0 &&
              ! is.nan(func_list[[sFuncName]]$argList[[paramName_pref]])) {
            paramName_ctrl = paramName_ctrl__default
          } else {
            I_found_value = FALSE
            for (argKey in listKeys_arg) {
              if (length(func_list[[sFuncName]]$argList[[argKey]]) != 0 &&
                  ! is.nan(func_list[[sFuncName]]$argList[[argKey]]) ) {
                I_set = TRUE
                paramName_ctrl = argKey
                break
              }
            }
          }

        }
      }
    }

    # Do we have a control param?
    if (is.null(paramName_ctrl)) {
      # Set to default
      paramName_ctrl = paramName_ctrl__default
    } else {
      # Validate paramName_ctrl is an accepted paramName
      func_list__valid = auction__valid_opt__func_list()
      param_list__valid = names(func_list__valid[[func_list[[sFuncName]]$funcName]])

      if (! paramName_ctrl %in% param_list__valid) {
        # Set to default
        paramName_ctrl = paramName_ctrl__default

        print(paste("auction__check_input__unobs_distrib__set_ctrl_arg",
                    " | Given control param is not valid for ",
                    sFuncName,
                    " - setting to default, ", paramName_ctrl__default,
                    sep=''))
      }
    }

    # Update 'func_list'
    #   set $argName_ctrl
    func_list[[sFuncName]]$argName_ctrl = paramName_ctrl
    #   set $argList to include ONLY $argName_ctrl
    if (is.list(func_list[[sFuncName]]$argList) && length(func_list[[sFuncName]]$argList) > 0) {
      # $argList is a list with entries
      #   is paramName_ctrl one of the entries?
      if (! paramName_ctrl %in% names(func_list[[sFuncName]]$argList) ||
          (is.nan(func_list[[sFuncName]]$argList[[paramName_ctrl]]) ||
           length(func_list[[sFuncName]]$argList[[paramName_ctrl]])==0 ) ) {
        func_list[[sFuncName]]$argList = list()
        func_list[[sFuncName]]$argList[[paramName_ctrl]] = NaN
      } else {
        # Remove other entires
        func_list[[sFuncName]]$argList = func_list[[sFuncName]]$argList[
          names(func_list[[sFuncName]]$argList) == paramName_ctrl ]
      }

    } else {
      # $argList is not a list or is an empty list
      func_list[[sFuncName]]$argList = list()
      func_list[[sFuncName]]$argList[[paramName_ctrl]] = NaN
    }

  }
  return( list( func_list=func_list, res=res ) )
}

auction__check_input__inital_guess__prep <- function(initial_guess, nParams_dat, res) {

  if ( !is.null(initial_guess) &&
       !is.vector(initial_guess) &&
       !is.list(initial_guess)
  ) {
    res['err_code'] = 2
    res['err_msg'] = "Invalid input for 'initial guess' | Must be NULL, a numeric vector, or a list"
    auction__gen_err_msg(res)
  }

  # Prepare 'initial_guess' with respect to:
  #   parameters unrelated to unobserved heterogeneity distribution
  if (! is.list(initial_guess) && ! is.character(initial_guess)) {
    # Generate default values for parameters unrelated to unobserved heterogeneity distribution
    initial_guess = auction__conv__init_pos__gen__without_unobs_distrib(nParams_dat=nParams_dat)

  } else if (is.vector(initial_guess)) {

    if (length(initial_guess) == nParams_dat) {
      # Initialize
      initial_guess__tmp = auction__conv__init_pos__gen__without_unobs_distrib(nParams_dat=nParams_dat)
      # Get corresponding indices
      listIdx = auction__x0_indices__without_unobs_distrib()
      # Fill from input 'initial_guess'
      initial_guess__tmp$pv_weibull_mu = initial_guess[listIdx$pv_weibull_mu]
      initial_guess__tmp$pv_weibull_a = initial_guess[listIdx$pv_weibull_a]
      initial_guess__tmp$x_terms = initial_guess[listIdx$x_terms__start:length(initial_guess)]

    } else if (length(names(func_list))==1 &&
               length(initial_guess)-1 == nParams_dat) {
      # Initialize
      initial_guess__tmp = auction__conv__init_pos__gen(nParams_dat=nParams_dat)
      # Get corresponding indices
      listIdx = auction__x0_indices()
      # Fill from input 'initial_guess'
      initial_guess__tmp$pv_weibull_mu = initial_guess[listIdx$pv_weibull_mu]
      initial_guess__tmp$pv_weibull_a = initial_guess[listIdx$pv_weibull_a]
      initial_guess__tmp$x_terms = initial_guess[listIdx$x_terms__start:length(initial_guess)]
      initial_guess__tmp$unobs_dist_param = initial_guess[listIdx$unobs_dist_param]

    } else {
      res['err_code'] = 2
      res['err_msg'] = "Not ready for vector input for 'initial guess'"
      auction__gen_err_msg(res)
    }
    # Done
    initial_guess = initial_guess__tmp
  }


  if (is.list(initial_guess)) {
    # Validate keys
    #   What keys are requied?
    listReq = c('pv_weibull_mu', 'pv_weibull_a', 'x_terms')
    #   What keys do we have?
    listKeys = names(initial_guess)
    #   Are there any missing keys?
    listMissing = listReq[! listReq %in% listKeys]
    if (length(listMissing) > 0) {
      res['err_code'] = 2
      res['err_msg'] = paste("Invalid input for 'initial_guess' | Missing required entries:",
                             paste(listMissing, collapse=', '))
      auction__gen_err_msg(res)
    }
    #   save to 'initial_guess__base'
    initial_guess__base = initial_guess[ listKeys %in% listReq ]
    #   Remove excess keys
    initial_guess__unobs = initial_guess[! listKeys %in% listReq ]

    if (! is.list(initial_guess__unobs) || length(initial_guess__unobs)==0) {
      initial_guess__unobs = list()
    }
  }
  return(list(
    base = initial_guess__base,
    unobs = initial_guess__unobs ))
}

auction__unobs_distrib__set_init_guess <- function(func_list, initial_guess) {
  # print("TEST - auction__unobs_distrib__set_init_guess")
  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)

    if (! is.nan(func_list[[sFuncName]]$argList[[func_list[[sFuncName]]$argName_ctrl]]) &&
        is.numeric(func_list[[sFuncName]]$argList[[func_list[[sFuncName]]$argName_ctrl]]) &&
        length(func_list[[sFuncName]]$argList[[func_list[[sFuncName]]$argName_ctrl]]) > 0) {
      func_list[[sFuncName]]$initial_guess = func_list[[sFuncName]]$argList[[func_list[[sFuncName]]$argName_ctrl]]
    } else {
      # Check to see if there is a name match within the variable 'initial_guess'
      I_match = FALSE
      if (! is.null(initial_guess) && ! is.list(initial_guess) && length(initial_guess) > 0) {
        listKeys = names(initial_guess)
        # Check if the entry within $funcName matches
        # then check if sFuncName matches
        if (func_list[[sFuncName]]$funcName %in% listKeys) {
          func_list[[sFuncName]]$initial_guess = initial_guess[[
            listKeys[ listKeys %in% func_list[[sFuncName]]$funcName ]
            ]]
          I_match = TRUE
        } else if (sFuncName %in% listKeys) {
          func_list[[sFuncName]]$initial_guess = initial_guess[[
            listKeys_initguess[ listKeys_initguess %in% sFuncName ]
            ]]
          I_match = TRUE
        }

      }
      # If no match found, get default
      if (! I_match) {
        func_list[[sFuncName]]$initial_guess =
          auction__conv__init_pos__gen__unobs_distrib(funcName = func_list[[sFuncName]]$funcName,
                                                      paramName = func_list[[sFuncName]]$argName_ctrl )
      }
    }


  }
  # print("TEST - auction__unobs_distrib__set_init_guess")
  return(func_list)
}
auction__unobs_distrib__set_step_size <- function(func_list) {
  for (funcName in names(func_list)) {
    sFuncName = as.character(funcName)

    if ( is.nan(func_list[[sFuncName]]$step_size) ||
         ! is.numeric(func_list[[sFuncName]]$step_size) ||
         length(func_list[[sFuncName]]$step_size)==0 ) {
      # Set value

      print(paste("Setting step-size for '", sFuncName, "' to 1", sep=''))

      func_list[[sFuncName]]$step_size = 1

    }
  }
  return(func_list)
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

auction__conv__init_pos__gen__unobs_distrib <- function(funcName, paramName = NULL) {
  if (funcName == 'dgamma') {
    if (! is.null(paramName) && (paramName=='rate' || paramName=='shape')) {
      if (paramName=='rate') {
        unobs_dist_param = 1
      } else {
        unobs_dist_param = 1
      }
    } else {
      unobs_dist_param = 1
    }

  } else if (funcName == 'dlnorm') {
    if (! is.null(paramName) && (paramName=='meanlog' || paramName=='sdlog')) {
      if (paramName=='meanlog') {
        unobs_dist_param = 0.5
      } else {
        unobs_dist_param = 0.5
      }
    } else {
      unobs_dist_param = 1
    }

  } else if (funcName == 'dweibull') {
    if (! is.null(paramName) && (paramName=='scale' || paramName=='shape')) {
      if (paramName=='scale') {
        unobs_dist_param = 1
      } else {
        unobs_dist_param = 1
      }
    } else {
      unobs_dist_param = 1
    }

  } else {
    unobs_dist_param = 1
  }
  return(unobs_dist_param)
}

auction__conv__init_pos__gen__without_unobs_distrib <- function(nParams_dat = NULL, nParams_X = 2) {
  # Goal: Generate x0
  #           x0 is the initial guess or initial position for
  #           the numerical solver
  #             -> optim()
  #
  #       Do not include parameter for Unobserved Heterogenity

  # Currently, fix initial guess
  pv_weibull_mu = 8
  pv_weibull_a = 2

  if (is.null(nParams_dat)) {
    h_x = rep(0.5, nParams_X)
  } else {
    # dat = "price" | "num of bids" | "x1" | "x2" | ...
    #   # X columns = # Columns - 2
    h_x = rep(0.5, (nParams_dat-2) )
  }
  return( list(
    pv_weibull_mu = pv_weibull_mu,
    pv_weibull_a = pv_weibull_a,
    x_terms = h_x
    ) )
}
auction__conv__init_pos__gen <- function(nParams_dat) {
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
  h_x = rep(0.5, (nParams_dat-2) )

  return( list(
    pv_weibull_mu = pv_weibull_mu,
    pv_weibull_a = pv_weibull_a,
    unobs_dist_param = unobs_dist_param,
    x_terms = h_x
    ) )
}

auction__conv__get_params <- function(func__unobs_distrib, initial_guess, max_iterations = 2000) {
  # Define convergence criteria
  #  Max iterations
  #  Step-size of parameters that may be adjusted
  # Define initial guess for convergence

  # Max iterations = maxit
  if (!is.nan(max_iterations) && is.numeric(max_iterations) && length(max_iterations) == 1) {
    maxit = max_iterations
  } else {
    maxit = 2000
  }

  # initial guess & step size
  #   Is parameter for unobserved heterogenity being held constant?
  I_unobs__const = (! is.nan(func__unobs_distrib$argList[[func__unobs_distrib$argName_ctrl]]) &&
                      is.numeric(func__unobs_distrib$argList[[func__unobs_distrib$argName_ctrl]]) &&
                      length(func__unobs_distrib$argList[[func__unobs_distrib$argName_ctrl]]) != 0 )
  #   Continue
  nParam = auction__get_num_columns__dat(initial_guess)
  nParam_X = length(initial_guess$x_terms)
  if  (! I_unobs__const) {
    nParam = nParam+1

    listIdx = auction__x0_indices()
    listStep = auction__x0_stepsizes()
  } else {
    listIdx = auction__x0_indices__without_unobs_distrib()
    listStep = auction__x0_stepsizes__without_unobs_distrib()
  }
  #   Initialize vectirs
  x0_init = rep(-1, nParam)
  x0_step = rep(-1, nParam)

  # Fill in base parameters (not unobserved heterogeneity distribution)
  x0_init[listIdx$pv_weibull_mu] = initial_guess$pv_weibull_mu
  x0_init[listIdx$pv_weibull_a] = initial_guess$pv_weibull_a
  x0_init[listIdx$x_terms__start:(listIdx$x_terms__start + nParam_X - 1)] = initial_guess$x_terms

  x0_step[listIdx$pv_weibull_mu] = listStep$pv_weibull_mu
  x0_step[listIdx$pv_weibull_a] = listStep$pv_weibull_a
  x0_step[listIdx$x_terms__start:(listIdx$x_terms__start + nParam_X - 1)] =
    c(listStep$x_terms__start, rep(listStep$x_terms__other, nParam_X-1))


  # If parameter for unobserved heterogeneity not held constant, add
  if (! I_unobs__const) {
    x0_init[listIdx$unobs_dist_param] = func__unobs_distrib$initial_guess
    x0_step[listIdx$unobs_dist_param] = func__unobs_distrib$step_size
    # x0_step[listIdx$unobs_dist_param] = listStep$unobs_dist_param
  }

  return( list( I_unobs__const = I_unobs__const,
                initial_guess = x0_init,
                conv_ctrl = list(maxit = max_iterations,
                                 parscale = x0_step )
                ) )
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
      shape=list(type='numeric',req=FALSE), # Shape is required
      scale=list(type='numeric',req=FALSE)
    )
  ))
}

auction__unobs_dist__ctrl_param__default <- function(funcName){
  if (funcName == 'dgamma') {
    paramName = 'rate'
  } else if (funcName == 'dlnorm') {
    paramName = 'sdlog'
  } else if (funcName == 'dweibull') {
    paramName = 'shape'
  } else {
    paramName = NULL
  }
}

auction__unobs_dist__exp_val_1 <- function(func__prob_distrib, paramVal) {
  # Goal: Update parameters within 'func__prob_distrib' to ensure
  #       E(X)=1 for the probability distribution representing
  #       unobserved heterogeneity

  paramName = func__prob_distrib$argName_ctrl
  if (func__prob_distrib$I_unobs__const) {
    paramVal = func__prob_distrib$argList[[paramName]]
  }

  # Ready to calculate value for other parameter
  if (func__prob_distrib$funcName == 'dgamma') {
    # Calculate
    if (paramName == 'rate') {
      list_paramVal = auction__unobs_dist__exp_val_1__gamma(param_rate=paramVal)
    } else {
      list_paramVal = auction__unobs_dist__exp_val_1__gamma(param_shape=paramVal)
    }
    # Store values with $argList
    if (is.null( names(func__prob_distrib$argList) )) {
      func__prob_distrib$argList = list(
        rate=list_paramVal$param_rate,
        shape=list_paramVal$param_shape
      )
    } else {
      func__prob_distrib$argList$rate = list_paramVal$param_rate
      func__prob_distrib$argList$shape = list_paramVal$param_shape
    }

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
    # Calculate
    if (paramName == 'shape') {
      list_paramVal = auction__unobs_dist__exp_val_1__weibull(param_shape=paramVal)
    } else {
      list_paramVal = auction__unobs_dist__exp_val_1__weibull(param_scale=paramVal)
    }
    # Store values with $argList
    if (is.null( names(func__prob_distrib$argList) )) {
      func__prob_distrib$argList = list(
        scale=list_paramVal$param_scale,
        shape=list_paramVal$param_shape
      )
    } else {
      func__prob_distrib$argList$scale = list_paramVal$param_scale
      func__prob_distrib$argList$shape = list_paramVal$param_shape
    }

  }
  return(func__prob_distrib)
}
auction__unobs_dist__std <- function(func__prob_distrib) {
  # Get standard deviation of unobserved heterogeneity
  if (func__prob_distrib$funcName == 'dgamma') {
    unobs_std = auction__unobs_dist__std__gamma(
      param_rate  = func__prob_distrib$argList$rate,
      param_shape = func__prob_distrib$argList$shape
    )
  } else if (func__prob_distrib$funcName == 'dlnorm') {
    unobs_std = auction__unobs_dist__std__lognorm(
      param_meanlog = func__prob_distrib$argList$meanlog,
      param_sdlog   = func__prob_distrib$argList$sdlog
    )
  } else if (func__prob_distrib$funcName == 'dweibull') {
    unobs_std = auction__unobs_dist__std__weibull(
      param_scale = func__prob_distrib$argList$scale,
      param_shape = func__prob_distrib$argList$shape
    )
  } else {
    print(paste("Unknown distribution '", func__prob_distrib$funcName,"', returning NaN for standard deviation", sep=''))
    unobs_std = NaN
  }
  return(unobs_std)
}

auction__unobs_dist__exp_val_1__lognorm <- function(param_meanlog = NULL, param_sdlog = NULL) {
  # Given either meanlog or sdlog, set the other such that E(x) = 1
  #   E(x) = exp(mu) * exp(1/2 * s ^ 2)

  if (is.null(param_meanlog) || is.nan(param_meanlog)) {
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

  if (! is.null(param_shape) || is.nan(param_shape)) {
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

  if (is.null(param_rate) || is.nan(param_rate)) {
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


  # Get appropriate list of indices for x0
  #   Is the parameter for unobserved heterogeneity held constant?
  if (func__prob_distrib$I_unobs__const) {
    listIdx = auction__x0_indices__without_unobs_distrib()
    unobs_paramVal = func__prob_distrib$argList[[func__prob_distrib$argName_ctrl]]
  } else {
    listIdx = auction__x0_indices()
    unobs_paramVal = x0[listIdx$unobs_dist_param]
  }


  h = x0[listIdx$x_terms__start:(listIdx$x_terms__start+dim(dat_X)[2]-1)]
  v__h = exp( colSums( h * t(dat_X) ) )

  if (! func__prob_distrib$I_unobs__const && x0[listIdx$unobs_dist_param] <= 0.1) {
    return(-Inf) # Check that these hold at estimated values
  } else if ( sum (x0[listIdx$pv_weibull_mu] <= 0 ) > 0 ) {
    return(-Inf)
  } else if ( sum( x0[listIdx$pv_weibull_a] <= 0.01 ) > 0) {
    return(-Inf)
  } else {
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/x0[listIdx$pv_weibull_a])
    v__w = dat_price / v__h

    dat = cbind(v__w, dat_num, x0[listIdx$pv_weibull_mu], x0[listIdx$pv_weibull_a], v__gamma_1p1opa)

    # Set E(X) = 1 for UnObserved distribution
    func__prob_distrib = auction__unobs_dist__exp_val_1(
      func__prob_distrib=func__prob_distrib,
      paramVal=unobs_paramVal )

    # Run
    v__f_w = parApply(cl = cl,
                      X = dat,
                      MARGIN = 1,
                      FUN = f__funk__v3,
                      func__prob_distrib=func__prob_distrib
    )

    # Return output
    ###print(paste("v__f_w", -sum(log(v__f_w / v__h))))
    ### print(func__prob_distrib)
    v__f_y = v__f_w / v__h
    return(-sum(log(v__f_y)))
  }
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
vf__bid_function_fast__v3 = Vectorize(FUN = f__bid_function_fast__v3,vectorize.args = "price")


auction__x0_indices <- function() {
  return( list(
    pv_weibull_mu = 1,
    pv_weibull_a = 2,
    unobs_dist_param = 3,
    x_terms__start = 4
    ) )
}
auction__x0_indices__without_unobs_distrib <- function() {
  return( list(
    pv_weibull_mu = 1,
    pv_weibull_a = 2,
    x_terms__start = 3
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
auction__x0_stepsizes__without_unobs_distrib <- function() {
  return( list(
    pv_weibull_mu = 1,
    pv_weibull_a = 0.1,
    x_terms__start = 0.1,
    x_terms__other = 1
  ) )
}


### CODE START

# log norm distribution
listInputPDF = list(dlnorm=list()) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)

#listInputPDF = lidlnorm=list()) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = "dlnorm")

# Weibull distribution
listInputPDF = list(dweibull=list()) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)

# Gamma distribution
listInputPDF = list(dgamma=list()) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)

# what happens if none specified?
listInputPDF = NULL # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)

# or what ahppens if wrong parameters specified?
listInputPDF = list(dlnorm=list(argList=list(meanlog=2)), dfake=list()) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)

# Lognorm again with fixed distribution of unobserved
listInputPDF = list(dlnorm=list(argList=list(sdlog=2))) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)

# run with multiple distributions requested
listInputPDF = list(dlnorm=list(), dweibull=list()) # working
auction_v2(dat = auction__generate_data(10),
           winning_bid = 'price',
           number_of_bids = 'num',
           num_cores = 3,
           func_list__unobs_distrib = listInputPDF)




stop()































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
