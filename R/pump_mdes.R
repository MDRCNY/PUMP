



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
#' @param tol tolerance for target power
#'
#' @param max.steps how many steps allowed before terminating
#' @param max.tnum maximum cumulative number of samples
#' @param final.tnum number of samples for final draw
#' @param cl cluster object to use for parallel processing
#' @param updateProgress the callback function to update the progress bar (User
#'   does not have to input anything)
#' @param just.result.table TRUE means only return final answer, FALSE means
#'   return search path information.
#'
#' @importFrom stats qt
#' @return mdes results
#' @export
#'

pump_mdes <- function(
  design, MTP = NULL, M, J, K = 1, numZero = NULL,
  target.power, power.definition, tol,
  nbar, Tbar, alpha,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  rho = NULL, rho.matrix = NULL, omega.2 = 0, omega.3 = 0,
  B = 1000,
  max.steps = 20, max.tnum = 2000, start.tnum = 200, final.tnum = 4*max.tnum,
  cl = NULL, updateProgress = NULL, give.optimizer.warnings = FALSE,
  just.result.table = TRUE,
  verbose = FALSE
)
{
  if ( verbose ) {
    scat( "pump_mdes with %d max iterations per search, starting at %d iterations with final %d iterations (%d perms for WY if used)\n",
          max.tnum, start.tnum, final.tnum, B )
  }

  if ( missing( "target.power" ) ||  missing( "power.definition" ) || missing( "tol" ) ) {
    stop( "target.power, power.definition, or tol (tolerance) not supplied" )
  }

  # validate input parameters
  params.list <- list(
    MTP = MTP, numZero = numZero,
    M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho, rho.matrix = rho.matrix, B = B
  )
  ##
  params.list <- validate_inputs(design, params.list, mdes.call = TRUE )
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
  mdes.cols <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))

  # check if zero power, then return 0 MDES
  if(round(target.power, 2) <= 0)
  {
    message('Target power of 0 (or negative) requested')
    mdes.results <- data.frame(MTP, 0, 0)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return( make.pumpresult( mdes.results,
                             type = "mdes",
                              params.list = params.list) )
  }

  # check if max power, then return infinite MDES
  if(round(target.power, 2) >= 1)
  {
    message('Target power of 1 (or larger) requested')
    mdes.results <- data.frame(MTP, Inf, 1)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return( make.pumpresult( mdes.results,
                             type = "mdes",
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
  Q.m <- calc.Q.m(
    design = design, J = J, K = K, nbar = nbar, Tbar = Tbar,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3
  )
  t.df <- calc.df(
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
  mdes.raw  <- ifelse(target.power > 0.5,
                      Q.m * (crit.alpha + crit.beta),
                      Q.m * (crit.alpha - crit.beta))
  mdes.bf   <- ifelse(target.power > 0.5,
                      Q.m * (crit.alphaxM + crit.beta),
                      Q.m * (crit.alphaxM - crit.beta))



  # MDES is already calculated for individual power for raw and Bonferroni
  if ( pdef$indiv & MTP == "Bonferroni") {
    mdes.results <- data.frame(MTP, mdes.bf, target.power)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results, type = "mdes",
                              params.list = params.list ) )
  }

  if ( MTP == "None") {
    mdes.results <- data.frame(MTP, mdes.raw, target.power)
    colnames(mdes.results) <- mdes.cols
    return( make.pumpresult( mdes.results, type = "mdes",
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

  optim.out <-optimize_power(design, search.type = 'mdes', MTP,
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

  test.pts <- optim.out$test.pts


  mdes.results <- data.frame(
    MTP,
    test.pts$pt[nrow(test.pts)],
    test.pts$power[nrow(test.pts)]
  )
  colnames(mdes.results) <- mdes.cols

  return( make.pumpresult( mdes.results, type = "mdes",
                           tries = test.pts,
                           params.list = params.list,
                           just.result.table = just.result.table  ) )
}




