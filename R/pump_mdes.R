#' MDES (minimum detectable effect size) function
#'
#' The minimum detectable effect size function calculates the most feasible
#' minimum detectable effect size for a given MTP, power and power definition.
#' The goal is to find the MDES value that satisfies the tolerance set in the
#' parameter in the power value.
#'
#' @inheritParams pump_power
#'
#' @param target.power Target power to arrive at
#' @param power.definition must be a valid power type outputted by power
#'   function, i.e. D1indiv, min1, etc.
#' @param tol tolerance for target power, defaults to 0.01 (1%).  This parameter
#'   controls when the search is done: when estimated power (checked with
#'   `final.tnum` iterations) is within `tol`, the search stops.
#'
#' @param max.steps how many steps allowed before terminating
#' @param start.tnum number of samples for first iteration of search algorithm
#' @param max.tnum maximum cumulative number of samples
#' @param final.tnum number of samples for final draw
#' @param cl cluster object to use for parallel processing
#' @param updateProgress the callback function to update the progress bar (User
#'   does not have to input anything)
#' @param give.optimizer.warnings whether to return verbose optimizer warnings
#'
#' @importFrom stats qt
#' @return mdes results
#' @export
#'

pump_mdes <- function(
  design, MTP = NULL, M, nbar, J, K = 1,
  Tbar, alpha = 0.05,
  target.power, power.definition, tol = 0.01,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  omega.2 = 0, omega.3 = 0,
  rho = NULL, rho.matrix = NULL,
  B = 1000,
  max.steps = 20, max.tnum = 2000, start.tnum = 200, final.tnum = 4*max.tnum,
  cl = NULL, updateProgress = NULL, give.optimizer.warnings = FALSE,
  verbose = FALSE
)
{

  # Call self for each element on MTP list.

  # NOTE: This is not well defined because do we store search history or what
  # when we have multiple calls to different MTPs?

  # if ( length( MTP ) > 1 ) {
  #   if ( verbose ) {
  #     scat( "Multiple MTPs leading to %d calls\n", length(MTP) )
  #   }
  #   des = purrr::map( MTP,
  #                     pump_mdes, design = design,
  #                     target.power = target.power, power.definition = power.definition, tol = tol,
  #                     M = M, J = J, K = K, nbar = nbar,
  #                     Tbar = Tbar, alpha = alpha,
  #                     numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
  #                     R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
  #                     ICC.2 = ICC.2, ICC.3 = ICC.3,
  #                     omega.2 = omega.2, omega.3 = omega.3,
  #                     rho = rho, rho.matrix = rho.matrix,
  #                     B = B,
  #                     max.steps = max.steps, max.tnum = max.tnum, start.tnum = start.tnum, final.tnum = final.tnum,
  #                     cl = cl, updateProgress = updateProgress, give.optimizer.warnings = give.optimizer.warnings,
  #                     verbose = verbose )
  #
  #   plist = attr( des[[1]], "params.list" )
  #   plist$MTP = MTP
  #     ftable = des[[1]]
  #     for ( i in 2:length(des) ) {
  #       ftable = dplyr::bind_rows( ftable, des[[i]] )
  #     }
  #
  #   return( make.pumpresult( ftable, "mdes",
  #                            params.list = plist,
  #                            design = design,
  #                            multiple_MTP = TRUE ) )
  #
  #   #des = map( des, ~ .x[nrow(.x),] ) %>%
  #   #  dplyr::bind_rows()
  #   #return( des )
  # }



  if ( verbose ) {
    scat( "pump_mdes with %d max iterations per search, starting at %d iterations with final %d iterations (%d perms for WY if used)\n",
          max.tnum, start.tnum, final.tnum, B )
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
    MTP = MTP,
    M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho, rho.matrix = rho.matrix, B = B
  )
  ##
  params.list <- validate_inputs(design, params.list, mdes.call = TRUE)
  ##
  MTP <- params.list$MTP
  MDES <- params.list$MDES
  M <- params.list$M; J <- params.list$J; K <- params.list$K
  nbar <- params.list$nbar; Tbar <- params.list$Tbar; alpha <- params.list$alpha
  numCovar.1 <- params.list$numCovar.1; numCovar.2 <- params.list$numCovar.2
  numCovar.3 <- params.list$numCovar.3
  R2.1 <- params.list$R2.1; R2.2 <- params.list$R2.2; R2.3 <- params.list$R2.3
  ICC.2 <- params.list$ICC.2; ICC.3 <- params.list$ICC.3
  omega.2 <- params.list$omega.2; omega.3 <- params.list$omega.3
  rho <- params.list$rho; rho.matrix <- params.list$rho.matrix
  B <- params.list$B

  # extract power definition
  pdef <- parse_power_definition( power.definition, M )

  # validate MTP
  if(MTP == 'None' & !pdef$indiv )
  {
    stop('For minimum or complete power, you must provide a MTP.')
  }

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
                             design = design,
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
                             design = design,
                             power.params.list = pow_params,
                             params.list = params.list) )
  }

  if ( verbose ) {
    message(paste("Estimating MDES for", MTP, "for target", power.definition,
                  "power of", round(target.power, 4)))
  }

  if (MTP == "WY-SD" && B < 1000){
    warning(paste("For the step-down Westfall-Young procedure,
                  it is recommended that sample (B) be at least 1000. Current B:", B))
  }

  # Compute Q.m and df
  Q.m <- calc_Q.m(
    design = design, J = J, K = K, nbar = nbar, Tbar = Tbar,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3
  )
  t.df <- calc_df(
    design = design, J = J, K = K, nbar = nbar,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3
  )

  # For raw and BF, compute critical values
  crit.alpha <- qt(p = (1-alpha/2), df = t.df)
  crit.alphaxM <- qt(p = (1-alpha/(2*M)), df = t.df)

  # Compute raw and BF MDES for individual power
  crit.beta <- ifelse(target.power > 0.5,
                      qt(target.power, df = t.df),
                      qt(1 - target.power, df = t.df))
  if(target.power > 0.5)
  {
    mdes.raw.list <- Q.m * (crit.alpha + crit.beta)
    mdes.bf.list  <- Q.m * (crit.alphaxM + crit.beta)
  } else
  {
    mdes.raw.list <- Q.m * (crit.alpha - crit.beta)
    mdes.bf.list  <- Q.m * (crit.alphaxM - crit.beta)
  }

  mdes.raw <- min(mdes.raw.list)
  mdes.bf <- max(mdes.bf.list)

  # MDES is already calculated for individual power for raw and Bonferroni
  if ( pdef$indiv & MTP == "Bonferroni") {
    mdes.results <- data.frame(MTP, mdes.bf, target.power)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results, type = "mdes",
                             design = design,
                             power.params.list = pow_params,
                             params.list = params.list ) )
  }

  if ( MTP == "None") {
    mdes.results <- data.frame(MTP, mdes.raw, target.power)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results, type = "mdes",
                             design = design,
                             power.params.list = pow_params,
                             params.list = params.list ) )
  }

  # MDES will be between raw and bonferroni for many power types
  mdes.low <- mdes.raw
  mdes.high <- mdes.bf

  # adjust bounds to capture needed range for minimum or complete power.
  # bounds note: complete power is a special case of minimum power
  if(pdef$min)
  {
    # complete power will have a higher upper bound
    # must detect all individual outcomes
    target.indiv.power <- target.power^(1/M)
    crit.beta <- ifelse(target.indiv.power > 0.5,
                        qt(target.indiv.power, df = t.df),
                        qt(1 - target.indiv.power, df = t.df))
    mdes.high   <- ifelse(target.indiv.power > 0.5,
                          Q.m * (crit.alphaxM + crit.beta),
                          Q.m * (crit.alphaxM - crit.beta))


    # min1 power will have a lower lower bound
    # must detect at least one individual outcome
    min.target.indiv.power <- 1 - (1 - target.power)^(1/M)
    crit.beta <- ifelse(min.target.indiv.power > 0.5,
                        qt(min.target.indiv.power, df = t.df),
                        qt(1 - min.target.indiv.power, df = t.df))
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

  test.pts <- optimize_power(design, search.type = 'mdes', MTP,
                             target.power, power.definition, tol,
                             start.tnum,
                             start.low = mdes.low, start.high = mdes.high,
                             MDES = NULL, J = J, K = K, nbar = nbar,
                             M = M, Tbar = Tbar, alpha = alpha,
                             numCovar.1 = numCovar.1,
                             numCovar.2 = numCovar.2,
                             numCovar.3 = numCovar.3,
                             R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                             ICC.2 = ICC.2, ICC.3 = ICC.3,
                             rho = rho, omega.2 = omega.2, omega.3 = omega.3,
                             B = B, cl = cl,
                             max.steps = max.steps, max.tnum = max.tnum,
                             final.tnum = final.tnum, give.warnings = give.optimizer.warnings)


  mdes.results <- data.frame(
    MTP,
    test.pts$pt[nrow(test.pts)],
    test.pts$power[nrow(test.pts)]
  )
  colnames(mdes.results) <- mdes.cols

  if(!is.na(mdes.results$`Adjusted.MDES`) && test.pts$dx[[nrow(test.pts)]] < 0.001 )
  {
    msg <- "Note: this function returns one possible value of MDES, but other (smaller values) may also be valid.\n"
    msg <- paste(msg, "Please refer to sample size vignette for interpretation.\n")
    message(msg)
    flat <- TRUE
  } else
  {
    flat <- FALSE
  }

  return( make.pumpresult( mdes.results, type = "mdes",
                           tries = test.pts,
                           flat = flat,
                           design = design,
                           params.list = params.list,
                           power.params.list = pow_params ) )
}




