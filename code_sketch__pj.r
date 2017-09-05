

auction <- function(
  input_data=NULL, initial_guess=NULL, generate_data=FALSE,
  input_pdf=NULL,
  num_cores=1
  ) {
  # INPUTS
  #   input_data
  #
  #   initial_guess [list]
  #
  #   generate_data [boolean]
  #
  #   input_pdf [string or list of strings]
  #     define the probability distribution function(s) to use
  #     input may be a string or list of strings
  #
  #   num_cores [int]
  #     Number of parallel workers to spawn.
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
  res = list(
    result=-1,
    err_code=0,
    err_msg=""
  )

  #   Validate environment
  #     Load required libraries
  res = auction__load_packages(res)

  # Check inputs
  #   Check parameter data types for non-null default value
  res = auction__check_input__generate_data(generate_data, res)
  res = auction__check_input__num_cores(num_cores, res)

  if ( res['err_code'] != 0 ) {
    return(res)
  }


  #   Check cores
  num_cores_detected = detectCores(res)
  if ( num_cores > num_cores_detected ) {
    print(paste0("Warning: You have requested ", num_cores,
                 " cores but only have ", num_cores_detected,
                 " cores available"))
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

  #   Check PDFs specified

  if ( res['err_code'] != 0 ) {
    return(res)
  }


  # Prepare data
  #   Part 1) Inspect data
  I_sampleData = FALSE
  if (generate_data) {
    I_sampleData = TRUE
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
  if ( I_sampleData ) {
    x0 = generate__data(v__y, v__n, m__h_x)
  } else {
    x0 = initial_guess
  }


  # Set up parallelization
  cl = makeCluster(num_cores)
  clusterExport(cl,varlist=c("vf__bid_function_fast",
                             "vf__w_integrand_z_fast",
                             "f__funk"))


  #   Get convergence conditions
  optim_control = f_conv_conditions(x0)


  # Run
  run_result = optim(par=x0, fn=f__ll_parallel, control=optim_control,
                 y=v__y, n=v__n, h_x=m__h_x, cl=cl)


  # Release resources
  stopCluster(cl)


  # Inspect result



  # Return result
  return(res)
}


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



generate__data <- function() {
  # For testing purposes, we will generate sample data

  set.seed(301)
  # data = # Generate some data
  # y, n, x1, x2: positive
  # n: discrete and > 1
  # y is some function of n, x1, x2

  obs = 200
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
vf__bid_function_fast = Vectorize(FUN = f.bid_function_fast,vectorize.args = "cost")

vf__w_integrand_z_fast = function(z, w_bid, num_bids, mu, alpha, gamma_1p1oa, param.u){

  b_z = vf__bid_function_fast(cost=z, num_bids=num_bids, mu=mu, alpha=alpha, gamma_1p1oa)
  u_z = w_bid/b_z

  vals = num_bids*alpha*(gamma_1p1oa/mu)^alpha*z^(alpha-1)*
    exp(-num_bids*(gamma_1p1oa/mu*z)^alpha)*
    1/b_z*
    dlnorm(u_z, meanlog=(-param.u^2*1/2), sdlog = param.u) # Note: can swap for different distributions

  vals[(gamma_1p1oa/mu)^alpha == Inf] = 0
  vals[exp(-num_bids*(gamma_1p1oa/mu*z)^alpha) == 0] = 0
  return(vals)

}



f__funk = function(data_vec, param_u){

  val = integrate(vf__w_integrand_z_fast, w_bid=data_vec[1],
                  num_bids=data_vec[2], mu=data_vec[3], alpha=data_vec[4],
                  gamma_1p1oa=data_vec[5], param_u=param_u, lower=0, upper=Inf,
                  abs.tol = 1e-10)
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
  else if ( sum( v.alpha <= 0.01 ) > 0)
    return(-Inf)
  else
    # Y Component
    v__gamma_1p1opa = gamma(1 + 1/v__alpha)
    v__w = v__y / v__h
    dat = cbind(v__w, v__n, v__mu, v__alpha, v__gamma_1p1opa)

    v__f_w = parApply(cl = cl, X = dat, MARGIN = 1, FUN = f__funk, param_u = u)
    v__f_y = v__f_w / v__h
    return(-sum(log(v.f_y)))
}






###########################################################################
# Run code
###########################################################################
res = auction(generate_data=TRUE, num_cores=20)
print(res)


