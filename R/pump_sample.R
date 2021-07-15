
# Code for the pump_sample method




#' This function calculates needed nbar to achieve a given power (as represented by
#' a difference in t-statistic) for all implemented designs
#'
#' @inheritParams calc.J
#'
#' @return nbar, the number of individuals needed, or NA if not possible given design

calc.nbar <- function(design, MT = 2.8, MDES, J, K = NULL, Tbar, R2.1,
                      R2.2, ICC.2, omega.2,
                      R2.3 = NULL, ICC.3 = NULL, omega.3 = NULL ) {

  if(design %in% c('blocked_i1_2c', 'blocked_i1_2f'))
  {
    nbar <- (MT/MDES)^2 * ( (1-ICC.2) * (1 - R2.1) / (Tbar * (1 - Tbar) * J) )
  } else if (design == 'blocked_i1_2r')
  {
    numr = (1-ICC.2)*(1-R2.1)
    denom = J * ((MDES/MT)^2) - ICC.2 * omega.2
    nbar <- numr / (Tbar*(1-Tbar)*denom)
  } else if (design == 'blocked_i1_3r') {
    numr = (1 - ICC.2 - ICC.3) * (1-R2.1)
    denom = J*K*((MDES/MT)^2) - J*ICC.3*omega.3 - ICC.2*omega.2
    nbar <- numr / ( Tbar*(1-Tbar)*denom )
  } else if (design == 'simple_c2_2r')
  {
    numr = (1-ICC.2)*(1-R2.1)
    denom = Tbar * (1-Tbar) * J * ((MDES/MT)^2) - ICC.2 * (1-R2.2)
    nbar <- numr / denom
  } else
  {
    stop(paste('Design not implemented:', design))
  }
  nbar <- ifelse( is.na( nbar ) || nbar < 0, NA, nbar )
  return( nbar )
}



#' This function calculates needed J to achieve a given power (as represented by
#' a difference in t-statistic) for all implemented designs
#'
#' @inheritParams pump_power
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MT Number of approximate effect-size unit SEs (adjusted for degrees of
#'   freedom issues) that the MDES needs to be to achieve desired power.  E.g.,
#'   2.8 for normal theory.
#' @param MDES scalar, or vector of length M; the MDES values for each outcome
#'
#' @return J, the number of schools needed

calc.J <- function(design, MT = 2.8, MDES, nbar, Tbar, R2.1, R2.2, ICC.2, omega.2) {

  if(design %in% c('blocked_i1_2c', 'blocked_i1_2f'))
  {
    J <- (MT/MDES)^2 * ( (1 - R2.1) / (Tbar * (1 - Tbar) * nbar) )
  } else if (design == 'blocked_i1_2r')
  {
    J <- (MT/MDES)^2 * ( (ICC.2 * omega.2) +
                           ((1 - ICC.2)*(1 - R2.1)) / (Tbar * (1 - Tbar) * nbar) )
  } else if (design == 'simple_c2_2r')
  {
    J <- (MT/MDES)^2 * ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar)) +
                           ((1 - ICC.2)*(1 - R2.1)) / (Tbar * (1 - Tbar) * nbar) )
  } else
  {
    stop(paste('Design not implemented:', design))
  }
  return( J )
}

#' Calculates K, the number of districts
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MT multiplier
#' @param vector of length M; the MDES values for each outcome.
#'   Can provide single MDES value which will be repeated for the M outcomes.
#' @param J scalar; the number of schools
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects
#'   variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects
#'   variability
#'
#' @return K, the number of districts
calc.K <- function(design, MT, MDES, J, nbar, Tbar,
                   R2.1, R2.2, R2.3,
                   ICC.2, ICC.3,
                   omega.2, omega.3) {

  if(design == 'blocked_i1_3r')
  {
    K <- (MT/MDES)^2 * ( (ICC.3 * omega.3) +
                           (ICC.2 * omega.2) / J +
                           ((1 - ICC.2 - ICC.3) * (1 - R2.1))/(Tbar * (1 - Tbar) * J * nbar) )
  } else if (design == 'simple_c3_3r')
  {
    K <- (MT/MDES)^2 * ( (ICC.3 * (1 - R2.3)) / (Tbar * (1 - Tbar)) +
                           (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
                           ((1 - ICC.2 - ICC.3)*(1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
  } else if (design == 'blocked_c2_3f')
  {
    K <- (MT/MDES)^2 * ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
                           ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
  } else if (design == 'blocked_c2_3r')
  {
    K <- (MT/MDES)^2 * ( (ICC.3 * omega.3) +
                           (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
                           ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
  } else
  {
    stop(paste('Design not implemented:', design))
  }
  return(K)
}





#' Calculating Sample for Raw (Unadjusted)
#'
#' This is a Helper function for getting Sample Size when no adjustments has
#' been made to the test statistics.
#'
#' It requires iteration because we do not know the degrees of freedom, and so
#' we guess and then calculate sample size, and then recalculate df based on
#' sample size, until we converge.
#'
#' @inheritParams pump_power
#'
#' @param design a single RCT design (see list/naming convention)
#' @param typesample type of sample size to calculate: J, K, or nbar
#' @param target.power target power to arrive at
#' @param tol tolerance of how much our search can change the sample size before
#'   stopping.  (In number of units, so 0.1 is relatively precise.)
#' @param two.tailed whether to calculate two-tailed or one-tailed power
#' @param max.steps how many steps allowed before terminating
#'
#' @return raw sample returns
#' @export

pump_sample_raw <- function(
  design, MTP, typesample,
  MDES,
  nbar = NULL, J = NULL, K = NULL,
  target.power,
  Tbar, alpha, two.tailed = TRUE,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2 = NULL, ICC.3 = NULL,
  omega.2 = NULL, omega.3 = NULL,
  tol = 0.1, max.steps = 100
)
{
  if ( typesample=="nbar" ) {
    stopifnot( is.null( nbar ) )
    nbar = 1000
  } else if ( typesample == "J" ) {
    stopifnot( is.null( J ) )
    J = 1000
  } else if ( typesample == "K" ) {
    stopifnot( is.null( K ) )
    K = 1000
  }

  i <- 0
  # convergence
  conv <- FALSE

  while (i <= max.steps & conv == FALSE) {
    # checking which type of sample we are estimating
    df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
    if(df < 0 | is.infinite(df)) {
      stop('Landed on situation with impossible df')
    }

    # t statistics
    T1 <- ifelse(two.tailed == TRUE, abs(qt(alpha/2, df)), abs(qt(alpha, df)))
    T2 <- abs(qt(target.power, df))

    # number of SEs we need for MDES
    MT <- ifelse(target.power >= 0.5, T1 + T2, T1 - T2)

    if (typesample == "J") {
      J1 <- calc.J(
        design, MT = MT, MDES = MDES[1], nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], ICC.2 = ICC.2[1], omega.2 = omega.2
      )
      if ( is.na(J1) || abs(J1 - J) < tol) {
        conv <- TRUE
      }
      J <- J1
    } else if (typesample == "K") {
      K1 <- calc.K(
        design, MT = MT, MDES = MDES[1], J = J, nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
        ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      if ( is.na(K1) || abs(K1 - K) < tol) {
        conv <- TRUE
      }
      K <- K1
    } else if (typesample == "nbar") {
      nbar1 <- calc.nbar(
        design, MT = MT, MDES = MDES[1], J = J, K = K, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
        ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      if (is.na( nbar1 ) || abs(nbar1 - nbar) < tol) {
        conv <- TRUE
      }
      nbar <- nbar1
    }
    i <- i + 1
  }

  if (typesample == "J") {
    return(J)
  } else if (typesample == "K") {
    return(K)
  } else if (typesample == "nbar") {
    return(nbar)
  }
}


parse_power_definition = function( power.definition, M ) {
  powertype = list( min = FALSE,
                    complete = FALSE,
                    indiv = FALSE )

  if ( stringr::str_detect( power.definition, "min" ) ) {
    powertype$min = TRUE
    powertype$min_k = as.numeric( gsub( "min", "", power.definition ) )
    stopifnot( !is.numeric( powertype$min_k ) )
  } else if ( stringr::str_detect( power.definition, "complete" ) ) {
    powertype$min = TRUE
    powertype$complete = TRUE
    powertype$min_k = M
  } else if ( stringr::str_detect( power.definition, "indiv" ) ) {
    powertype$indiv = TRUE
    powertime$indiv_k = as.numeric( gsub( "indiv", "", power.definition ) )
    stopifnot( !is.numeric( powertype$indiv_k ) )
  }

  return( powertype )
}


#' Calculate sample size
#'
#' Note: These currently only work if MDES is the same for all outcomes.
#'
#' @inheritParams pump_power
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported
#'   options: Bonferroni, BH, Holm, WY-SS, WY-SD
#' @param typesample type of sample size to calculate: "J", "K", or "nbar".
#' @param MDES scalar, or vector of length M; the MDES values for each outcome.
#' @param target.power target power to arrive at
#' @param power.definition must be a valid power type output by power function,
#'   i.e. D1indiv, min1, etc.
#' @param tol tolerance
#' @param two.tailed whether to calculate two-tailed or one-tailed power
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
#' @return sample size results
#' @export

pump_sample <- function(
  design, MTP, typesample,
  MDES, M,
  nbar = NULL, J = NULL, K = NULL,
  target.power, power.definition, tol,
  alpha, two.tailed = TRUE,
  Tbar,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  rho,
  omega.2 = 0, omega.3 = 0,
  tnum = 10000, B = 1000,
  max.steps = 20, max.cum.tnum = 5000, start.tnum = 200, final.tnum = 10000,
  cl = NULL, updateProgress = NULL
)
{
  # Give prelim values for the validation of parameters process.
  if ( typesample=="nbar" ) {
    stopifnot( is.null( nbar ) )
    nbar = 1000
  } else if ( typesample == "J" ) {
    stopifnot( is.null( J ) )
    J = 1000
  } else if ( typesample == "K" ) {
    stopifnot( is.null( K ) )
    K = 1000
  }

  # validation
  if(length(MDES) > 1 & length(unique(MDES)) > 1)
  {
    stop('Procedure assumes MDES is the same for all outcomes.')
  }

  # validate input parameters
  params.list <- list(
    MDES = MDES, M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho
  )
  ##
  params.list <- pum:::validate_inputs(design, MTP, params.list,single_MDES=TRUE)
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

  # Delete parameter we are actually going to search over.
  if ( typesample=="nbar" ) {
    nbar = NULL
  } else if ( typesample == "J" ) {
    J = NULL
  } else if ( typesample == "K" ) {
    K = NULL
  }

  output.colnames <- c("MTP", "Sample type", "Sample size",
                       paste(power.definition, "power") )

  # check if zero power
  if(round(target.power, 2) == 0)
  {
    message('Target power of 0 requested')
    test.pts <- NULL
    ss.results <- data.frame(MTP, typesample, 0, 0)
    colnames(ss.results) <- output.colnames
    return(list(ss.results = ss.results, test.pts = test.pts))
  }

  # Checks on what we are estimating, sample size
  message(paste("Estimating sample size of type", typesample, "for",
                MTP, "for target",
                power.definition, "power of", round(target.power, 4)))

  # Progress Message for the Type of Sample we are estimating, the type of power
  # and the targeted power value
  if (is.function(updateProgress)) {
    msg <- (paste("Estimating", whichSS, "for target", power.definition,
                  "power of",round(power,4)))
    updateProgress(message = msg)
  }


  # Compute needed sample size for raw and BF SS for INDIVIDUAL POWER. We are
  # estimating (potential) bounds like we estimated MDES bounds.
  #
  # For now assuming only two tailed tests
  ss.raw <- pump_sample_raw(
    design = design, MTP=MTP, typesample=typesample,
    MDES=MDES, J=J, K=K,
    target.power=target.power,
    nbar=nbar, Tbar=Tbar, alpha=alpha, two.tailed=two.tailed,
    numCovar.1=numCovar.1, numCovar.2=numCovar.2, numCovar.3=numCovar.3,
    R2.1=R2.1, R2.2=R2.2, R2.3=R2.3, ICC.2=ICC.2, ICC.3=ICC.3,
    omega.2=omega.2, omega.3=omega.3 )

  # We are done if raw power is what we are looking for
  if (MTP == "rawp"){
    raw.ss <- data.frame(MTP, power.definition, ss.raw, typesample, target.power)
    colnames(raw.ss) <- output.colnames
    return(raw.ss)
  }

  # Now identify sample size for Bonferroni
  ss.BF <- pump_sample_raw(
    design = design, MTP=MTP, typesample=typesample,
    MDES=MDES, J=J, K=K,
    target.power=target.power,
    nbar=nbar, Tbar=Tbar,
    alpha=alpha / M, # change alpha for BF
    two.tailed=two.tailed,
    numCovar.1=numCovar.1, numCovar.2=numCovar.2, numCovar.3=numCovar.3,
    R2.1=R2.1, R2.2=R2.2, R2.3=R2.3, ICC.2=ICC.2, ICC.3=ICC.3,
    omega.2=omega.2, omega.3=omega.3 )

  # Done if Bonferroni is what we are looking for
  if (MTP == "Bonferroni") {
    ss.BF <- data.frame(MTP, power.definition, ss.BF, typesample, target.power)
    colnames(ss.BF) <- output.colnames
    return(ss.BF)
  }

  # Like the MDES calculation, the sample size would be between raw and Bonferroni.
  # There is no adjustment and there is very conservative adjustment
  ss.low <- ss.raw
  ss.high <- ss.BF

  pdef = parse_power_definition( power.definition, M )

  if ( pdef$min ) {
    # adjust bounds to capture needed range (assuming independence, so approximate)
    need_pow = 1 - (1 - target.power)^(1/M)
    ss.low <- pump_sample_raw(
      design = design, MTP=MTP, typesample=typesample,
      MDES=MDES, J=J, K=K,
      target.power=need_pow,
      nbar=nbar, Tbar=Tbar,
      alpha=alpha / M, # change alpha for BF
      two.tailed=two.tailed,
      numCovar.1=numCovar.1, numCovar.2=numCovar.2, numCovar.3=numCovar.3,
      R2.1=R2.1, R2.2=R2.2, R2.3=R2.3, ICC.2=ICC.2, ICC.3=ICC.3,
      omega.2=omega.2, omega.3=omega.3 )

    need_pow = (target.power^(1/M))
    ss.high <- pump_sample_raw(
      design = design, MTP=MTP, typesample=typesample,
      MDES=MDES, J=J, K=K,
      target.power=need_pow,
      nbar=nbar, Tbar=Tbar,
      alpha=alpha / M, # change alpha for BF
      two.tailed=two.tailed,
      numCovar.1=numCovar.1, numCovar.2=numCovar.2, numCovar.3=numCovar.3,
      R2.1=R2.1, R2.2=R2.2, R2.3=R2.3, ICC.2=ICC.2, ICC.3=ICC.3,
      omega.2=omega.2, omega.3=omega.3 )
  }



  # If we can't make it work with raw, then we can't make it work.
  if ( is.na( ss.low ) ) {
    ss <- data.frame(MTP, power.definition, NA, typesample, target.power, typesample)
    colnames(ss) <- output.colnames
    return(ss)
  }


  ss.low = round( ss.low )
  ss.high = round( ss.high )



  # sometimes we already know the answer!
  if(ss.low == ss.high)
  {
    test.pts <- NULL
    ss.results <- data.frame(MTP, typesample, 1, target.power)
    colnames(ss.results) <- output.colnames
    return(list(ss.results = ss.results, test.pts = test.pts))
  }

  # search in the grid from min to max.
  test.pts <- optimize_power(
    design = design, search.type = typesample,
    MTP, target.power, power.definition, tol,
    start.tnum, start.low = ss.low, start.high = ss.high,
    MDES = MDES,
    J = J, K = K, nbar = nbar,
    M = M, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
    rho = rho, omega.2 = omega.2, omega.3 = omega.3,
    B = B, cl = cl,
    max.steps = max.steps, max.cum.tnum = max.cum.tnum,
    final.tnum = final.tnum
  )

  # Assemble results
  ss.results <- data.frame(
    MTP,
    typesample,
    ifelse(is.na(test.pts$pt[nrow(test.pts)]),
           NA,   # failed to find solution
           ceiling(test.pts$pt[nrow(test.pts)])),  # round up to get nice sufficient sample size.
    test.pts$power[nrow(test.pts)]
  )
  colnames(ss.results) <- output.colnames

  return(list(ss.results = ss.results, test.pts = test.pts))
}


