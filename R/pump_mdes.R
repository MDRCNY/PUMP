#' @title Estimate the minimum detectable effect size (MDES) (core function)
#'
#' @description The user chooses the context (d_m), MTP,
#' power definition, and choices of all relevant design parameters.
#' 
#' The functions performs a search algorithm,
#' and returns the MDES value within the specified tolerance.
#' For a list of choices for specific parameters, see pump_info().
#' 
#' 
#' @seealso For more detailed information about this function 
#' and the user choices,
#' see the manuscript \url{https://arxiv.org/abs/2112.15273},
#' which includes a detailed Technical Appendix
#' including information about the designs and models
#' and parameters.
#'
#' @inheritParams pump_power
#'
#' @param target.power target power for search algorithm.
#' @param power.definition see pump_info() for 
#' possible power definitions.
#' @param tol tolerance for target power, defaults to 0.01 (1%).  
#' This parameter controls when the search is done: 
#' when estimated power (checked with `final.tnum` iterations) 
#' is within `tol`, the search stops.
#' @param max.steps how many steps allowed before terminating.
#' @param tnum max number of samples for first iteration 
#' of search algorithm.
#' @param start.tnum number of samples to start search 
#' (this will increase with each step).
#' @param final.tnum number of samples for final draw.
#' @param give.optimizer.warnings whether to return 
#' verbose optimizer warnings.
#'
#' @return a pumpresult object containing MDES results.
#' @export
#'
#' @examples 
#' mdes <-  pump_mdes(
#'   d_m = "d3.1_m3rr2rr",
#'   MTP = 'HO',
#'   power.definition = 'D1indiv',
#'   target.power = 0.6,
#'   J = 30,
#'   K = 15,
#'   nbar = 50,
#'   M = 3,
#'   Tbar = 0.5, alpha = 0.05, 
#'   two.tailed = FALSE,
#'   numCovar.1 = 1, numCovar.2 = 1,
#'   R2.1 = 0.1, R2.2 = 0.1,
#'   ICC.2 = 0.2, ICC.3 = 0.2,
#'   omega.2 = 0.1, omega.3 = 0.1, 
#'   rho = 0.5, tnum = 2000)

pump_mdes <- function(
  d_m, MTP = NULL, numZero = NULL, M, nbar, J, K = 1,
  Tbar, alpha = 0.05, two.tailed = TRUE,
  target.power, power.definition, tol = 0.01,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  omega.2 = 0, omega.3 = 0,
  rho = NULL, rho.matrix = NULL,
  B = 1000,
  max.steps = 20, 
  tnum = 1000, start.tnum = tnum / 10, final.tnum = 4*tnum,
  parallel.WY.cores = 1,
  updateProgress = NULL, give.optimizer.warnings = FALSE,
  verbose = FALSE
)
{
  if ( verbose ) {
    scat( "pump_mdes with %d max iterations per search, 
          starting at %d iterations with final %d iterations 
          (%d perms for WY if used)\n",
          start.tnum, tnum, final.tnum, B )
  }

  if ( missing( "target.power" ) ||  missing( "power.definition" )  ) {
    stop( "target.power or power.definition not supplied" )
  }
  if ( is.null( "tol" ) ) {
    stop( "Cannot have NULL tol (tolerance)" )
  }

  pow_params <- list( target.power = target.power,
                      power.definition = power.definition,
                      tol = tol )

  # validate input parameters
  params.list <- list(
    MTP = MTP, numZero = numZero, 
    M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha, two.tailed = two.tailed,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho, rho.matrix = rho.matrix, B = B,
    max.steps = max.steps, 
    start.tnum = start.tnum, tnum = tnum, final.tnum = final.tnum,
    power.definition = power.definition
  )
  ##
  params.list <- validate_inputs(
      d_m, params.list, mdes.call = TRUE, verbose = verbose 
  )
  ##
  MTP <- params.list$MTP
  MDES <- params.list$MDES; numZero <- params.list$numZero
  M <- params.list$M; J <- params.list$J; K <- params.list$K
  nbar <- params.list$nbar; Tbar <- params.list$Tbar
  alpha <- params.list$alpha; two.tailed <- params.list$two.tailed
  numCovar.1 <- params.list$numCovar.1; numCovar.2 <- params.list$numCovar.2
  numCovar.3 <- params.list$numCovar.3
  R2.1 <- params.list$R2.1; R2.2 <- params.list$R2.2; R2.3 <- params.list$R2.3
  ICC.2 <- params.list$ICC.2; ICC.3 <- params.list$ICC.3
  omega.2 <- params.list$omega.2; omega.3 <- params.list$omega.3
  rho <- params.list$rho; rho.matrix <- params.list$rho.matrix
  B <- params.list$B
  power.definition <- params.list$power.definition
  params.list <- params.list[names(params.list) != 'power.definition']

  # extract power definition
  pdef <- parse_power_definition( power.definition, M )

  # information that will be returned to the user
  mdes.cols <- c("MTP", "Adjusted.MDES", paste(power.definition, "power"))

  # check if zero power, then return 0 MDES
  if(round(target.power, 2) <= 0)
  {
    message('Target power of 0 (or negative) requested')
    mdes.results <- data.frame(MTP, 0, 0)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results,
                             type = "mdes",
                             d_m = d_m,
                             power.params.list = pow_params,
                             params.list = params.list) )
  }

  # check if max power, then return infinite MDES
  if(round(target.power, 2) >= 1)
  {
    message('Target power of 1 (or larger) requested')
    mdes.results <- data.frame(MTP, Inf, 1)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results,
                             type = "mdes",
                             d_m = d_m,
                             power.params.list = pow_params,
                             params.list = params.list) )
  }

  if ( verbose ) {
    message(paste("Estimating MDES for", MTP, "for target", power.definition,
                  "power of", round(target.power, 4)))
  }

  # Compute Q.m and df
  Q.m <- calc_SE(
    d_m = d_m, J = J, K = K, nbar = nbar, Tbar = Tbar,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3
  )
  t.df <- calc_df(
    d_m = d_m, J = J, K = K, nbar = nbar,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3
  )

  # For raw and BF, compute critical values
  if(two.tailed)
  {
    crit.alpha <- stats::qt(p = (1-alpha/2), df = t.df)
    crit.alphaxM <- stats::qt(p = (1-alpha/(2*M)), df = t.df)
  } else
  {
    crit.alpha <- stats::qt(p = (1-alpha), df = t.df)
    crit.alphaxM <- stats::qt(p = (1-alpha/M), df = t.df)
  }


  # Compute raw and BF MDES for individual power
  crit.beta <- ifelse(target.power > 0.5,
                      stats::qt(target.power, df = t.df),
                      stats::qt(1 - target.power, df = t.df))
  if(target.power > 0.5)
  {
    mdes.raw.list <- Q.m * (crit.alpha + crit.beta)
    mdes.bf.list  <- Q.m * (crit.alphaxM + crit.beta)
  } else
  {
    mdes.raw.list <- Q.m * (crit.alpha - crit.beta)
    mdes.bf.list  <- Q.m * (crit.alphaxM - crit.beta)
  }

  mdes.low <- min(mdes.raw.list)
  mdes.high <- max(mdes.bf.list)

  # MDES is already calculated for individual power for raw and Bonferroni
  if ( pdef$indiv & MTP == "BF") {
    mdes.results <- data.frame(MTP, mdes.high, target.power)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results, type = "mdes",
                             d_m = d_m,
                             power.params.list = pow_params,
                             params.list = params.list ) )
  }

  if ( MTP == "None") {
    mdes.results <- data.frame(MTP, mdes.low, target.power)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results, type = "mdes",
                             d_m = d_m,
                             power.params.list = pow_params,
                             params.list = params.list ) )
  }

  # adjust bounds to capture needed range for minimum or complete power.
  # bounds note: complete power is a special case of minimum power
  if(pdef$min)
  {
    # complete power will have a higher upper bound
    # must detect all individual outcomes
    target.indiv.power <- target.power^(1/M)
    crit.beta <- ifelse(target.indiv.power > 0.5,
                        stats::qt(target.indiv.power, df = t.df),
                        stats::qt(1 - target.indiv.power, df = t.df))
    mdes.high   <- ifelse(target.indiv.power > 0.5,
                          Q.m * (crit.alphaxM + crit.beta),
                          Q.m * (crit.alphaxM - crit.beta))


    # min1 power will have a lower lower bound
    # must detect at least one individual outcome
    min.target.indiv.power <- 1 - (1 - target.power)^(1/M)
    crit.beta <- ifelse(min.target.indiv.power > 0.5,
                        stats::qt(min.target.indiv.power, df = t.df),
                        stats::qt(1 - min.target.indiv.power, df = t.df))
    mdes.low  <- ifelse(min.target.indiv.power > 0.5,
                        Q.m * (crit.alpha + crit.beta),
                        Q.m * (crit.alpha - crit.beta))


  }

  # unlikely, but just in case
  if(mdes.high < 0)
  {
    mdes.high <- 0
  }
  # corner case
  if(mdes.low < 0)
  {
    mdes.low <- 0
  }

  test.pts <- optimize_power(d_m, search.type = 'mdes', MTP,
                             target.power, power.definition, tol,
                             start.low = mdes.low, start.high = mdes.high,
                             MDES = NULL, J = J, K = K, nbar = nbar,
                             M = M, numZero = numZero, Tbar = Tbar, 
                             alpha = alpha, 
                             two.tailed = two.tailed,
                             numCovar.1 = numCovar.1,
                             numCovar.2 = numCovar.2,
                             numCovar.3 = numCovar.3,
                             R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                             ICC.2 = ICC.2, ICC.3 = ICC.3,
                             rho = rho, omega.2 = omega.2, omega.3 = omega.3,
                             B = B, parallel.WY.cores = parallel.WY.cores,
                             max.steps = max.steps, 
                             tnum = tnum, start.tnum = start.tnum, 
                             final.tnum = final.tnum, 
                             give.warnings = give.optimizer.warnings)


  mdes.results <- data.frame(
    MTP,
    test.pts$pt[nrow(test.pts)],
    test.pts$power[nrow(test.pts)]
  )
  colnames(mdes.results) <- mdes.cols

  if(!is.na(mdes.results$`Adjusted.MDES`) && 
     test.pts$dx[[nrow(test.pts)]] < 0.001 )
  {
    msg <- "Note: this function returns one possible value of MDES, 
    but other (smaller values) may also be valid.\n"
    msg <- paste(msg, "Please refer to sample size vignette for 
                 interpretation.\n")
    message(msg)
    flat <- TRUE
  } else
  {
    flat <- FALSE
  }

  return( make.pumpresult( mdes.results, type = "mdes",
                           tries = test.pts,
                           flat = flat,
                           d_m = d_m,
                           params.list = params.list,
                           power.params.list = pow_params ) )
}




