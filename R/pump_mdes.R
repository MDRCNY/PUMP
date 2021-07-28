
#' MDES (minimum detectable effect size) function
#'
#' The minimum detectable effect size function calculates the most feasible
#' minimum detectable effect size for a given MTP, power and power definition.
#' The goal is to find the MDES value that satisfies the tolerance set in the
#' parameter in the power value.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported
#'   options: Bonferroni, BH, Holm, WY-SS, WY-SD
#' @param target.power Target power to arrive at
#' @param power.definition must be a valid power type outputted by power
#'   function, i.e. D1indiv, min1, etc.
#' @param tol tolerance for target power
#' @param M scalar; the number of hypothesis tests (outcomes)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the
#'   treatment
#' @param alpha scalar; the family wise error rate (FWER)
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not
#'   including block dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#' @param R2.1 scalar, or vector of length M; percent of variation explained by
#'   Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by
#'   Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by
#'   Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random
#'   effects variability
#' @param omega.3 scalar; ratio of district effect size variability to random
#'   effects variability
#' @param rho scalar; correlation between outcomes
#' @param tnum scalar; the number of test statistics (samples)
#' @param B scalar; the number of samples/permutations for Westfall-Young
#' @param max.steps how many steps allowed before terminating
#' @param max.cum.tnum maximum cumulative number of samples
#' @param final.tnum number of samples for final draw
#' @param cl cluster object to use for parallel processing
#' @param updateProgress the callback function to update the progress bar (User
#'   does not have to input anything)
#'
#' @importFrom stats qt
#' @return mdes results
#' @export
#'

pump_mdes <- function(
  design, MTP, M, J, K = 1,
  target.power, power.definition, tol,
  nbar, Tbar, alpha,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  rho, omega.2 = 0, omega.3 = 0,
  tnum = 10000, B = 1000,
  max.steps = 20, max.cum.tnum = 5000, start.tnum = 200, final.tnum = 10000,
  cl = NULL, updateProgress = NULL, give.optimizer.warnings = FALSE
)
{
  if ( missing( "target.power" ) ||  missing( "power.definition" ) || missing( "tol" ) ) {
    stop( "target.power, power.definition, or tol (tolerance) not supplied" )
  }

  # validate input parameters
  params.list <- list(
    M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho
  )
  ##
  params.list <- validate_inputs(design, MTP, params.list, mdes.call = TRUE )
  ##
  MDES <- params.list$MDES
  M <- params.list$M; J <- params.list$J; K <- params.list$K
  nbar <- params.list$nbar; Tbar <- params.list$Tbar; alpha <- params.list$alpha
  numCovar.1 <- params.list$numCovar.1; numCovar.2 <- params.list$numCovar.2
  numCovar.3 <- params.list$numCovar.3
  R2.1 <- params.list$R2.1; R2.2 <- params.list$R2.2; R2.3 <- params.list$R2.3
  ICC.2 <- params.list$ICC.2; ICC.3 <- params.list$ICC.3
  omega.2 <- params.list$omega.2; omega.3 <- params.list$omega.3
  rho <- params.list$rho

  # check if zero power, then return 0 MDES
  if(round(target.power, 2) <= 0)
  {
    message('Target power of 0 (or negative) requested')
    test.pts <- NULL
    mdes.results <- data.frame(MTP, 0, 0)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return(list(mdes.results = mdes.results, test.pts = test.pts))
  }

  # check if max power, then return infinite MDES
  if(round(target.power, 2) >= 1)
  {
    message('Target power of 1 (or larger) requested')
    test.pts <- NULL
    mdes.results <- data.frame(MTP, Inf, 1)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return(list(mdes.results = mdes.results, test.pts = test.pts))
  }

  message(paste("Estimating MDES for", MTP, "for target", power.definition,
                "power of", round(target.power, 4)))

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
  crit.alphaxM <- qt(p = (1-alpha/M/2), df = t.df)

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

  pdef <- parse_power_definition( power.definition, M )

  # MDES is alrady calculated for individual power for raw and Bonferroni
  if ( pdef$indiv ) {
    if (MTP == "rawp"){
      mdes.results <- data.frame(MTP, mdes.raw, target.power)
      colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
      return (list(mdes.results = mdes.results, tries = NULL))
    } else if (MTP == "Bonferroni"){
      mdes.results <- data.frame(MTP, mdes.bf, target.power)
      colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
      return(list(mdes.results = mdes.results, tries = NULL))
    }
  }

  # complete power
  if(pdef$complete)
  {
    # must detect all individual outcomes
    target.indiv.power <- target.power^(1/M)
    crit.beta <- ifelse(target.indiv.power > 0.5,
                        qt(target.indiv.power, df = t.df),
                        qt(1 - target.indiv.power, df = t.df))
    mdes.bf   <- ifelse(target.indiv.power > 0.5,
                        Q.m * (crit.alphaxM + crit.beta),
                        Q.m * (crit.alphaxM - crit.beta))
  }

  # min power
  if(pdef$min)
  {
    # min1 power is going to be a lower bound
    # must detect at least one individual outcome
    min.target.indiv.power <- 1 - (1 - target.power)^(1/M)
    crit.beta <- ifelse(min.target.indiv.power > 0.5,
                        qt(min.target.indiv.power, df = t.df),
                        qt(1 - min.target.indiv.power, df = t.df))
    mdes.raw  <- ifelse(min.target.indiv.power > 0.5,
                        Q.m * (crit.alpha + crit.beta),
                        Q.m * (crit.alpha - crit.beta))
  }

  # MDES will be between raw and bonferroni
  mdes.low <- mdes.raw
  mdes.high <- mdes.bf

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
                             max.steps = max.steps, max.cum.tnum = max.cum.tnum,
                             final.tnum = final.tnum, give.warnings = give.optimizer.warnings)

  mdes.results <- data.frame(MTP, test.pts$pt[nrow(test.pts)], test.pts$power[nrow(test.pts)])
  colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))

  return(list(mdes.results = mdes.results, test.pts = test.pts))

}



