
# Code for the pump_sample method




#' This function calculates needed nbar to achieve a given power (as represented by
#' a difference in t-statistic) for all implemented designs
#'
#' @inheritParams calc.J
#'
#' @return nbar, the number of individuals needed, or NA if not possible given design
#' @export

calc.nbar <- function(design, MT = 2.8, MDES, J, K = NULL, Tbar, R2.1,
                      R2.2, ICC.2, omega.2,
                      R2.3 = NULL, ICC.3 = NULL, omega.3 = NULL ) {

  if(design %in% c('d1.1_m2cc'))
  {
    numr <- (1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J
    nbar <- (MT/MDES)^2 * numr/denom
  } else if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
  {
    numr <- (1 - ICC.2) * (1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J
    nbar <- (MT/MDES)^2 * numr/denom
  } else if (design == 'd2.1_m2fr')
  {
    numr <- (1 - ICC.2)*(1 - R2.1)
    denom <- J * ((MDES/MT)^2) - ICC.2 * omega.2
    nbar <- numr / (Tbar*(1-Tbar)*denom)
  } else if (design == 'd3.1_m3rr2rr') {
    numr <- (1 - ICC.2 - ICC.3) * (1 - R2.1)
    denom <- J*K*((MDES/MT)^2) - J*ICC.3*omega.3 - ICC.2*omega.2
    nbar <- numr / ( Tbar*(1-Tbar)*denom )
  } else if (design == 'd2.2_m2rc')
  {
    numr <- (1 - ICC.2)*(1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J * ((MDES/MT)^2) - ICC.2 * (1 - R2.2)
    nbar <- numr / denom
  } else if (design == 'd3.3_m3rc2rc')
  {
    numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J * K * ((MDES/MT)^2) - J * ICC.3*(1 - R2.3)  - ICC.2 * (1 - R2.2)
    nbar <- numr / denom
  } else if (design == 'd3.2_m3ff2rc')
  {
    numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J * K * ((MDES/MT)^2) - ICC.2 * (1 - R2.2)
    nbar <- numr / denom
  } else if (design == 'd3.2_m3rr2rc')
  {
    numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J * ( K * ((MDES/MT)^2) - ICC.3 * omega.3 ) - ICC.2 * (1 - R2.2)
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
#' @export

calc.J <- function(
  design, MT = 2.8, MDES, K = NULL, nbar, Tbar,
  R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3
) {

  if(design %in% c('d1.1_m2cc'))
  {
    numr <- (1 - R2.1)
    denom <- (Tbar * (1 - Tbar) * nbar)
    J <- (MT/MDES)^2 * numr/denom
  } else if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
  {
    numr <- (1 - ICC.2) * (1 - R2.1)
    denom <- (Tbar * (1 - Tbar) * nbar)
    J <- (MT/MDES)^2 * numr/denom
  } else if (design == 'd2.1_m2fr')
  {
    numr <- (1 - ICC.2) * (1 - R2.1)
    denom <- (Tbar * (1 - Tbar) * nbar)
    J <- (MT/MDES)^2 * ( (ICC.2 * omega.2) + numr / denom)
  } else if (design == 'd3.1_m3rr2rr')
  {
    numr <- (1 - ICC.2 - ICC.3 ) * (1 - R2.1) + Tbar * (1 - Tbar) * nbar * ICC.2 * omega.2
    denom <- K * (MDES/MT)^2 - ICC.3 * omega.3
    J <- (1 / (Tbar * (1 - Tbar) * nbar)) * numr/denom
  } else if (design == 'd2.2_m2rc')
  {
    numr <- nbar * ICC.2 * (1 - R2.2) + (1 - ICC.2) * (1 - R2.1)
    denom <- Tbar * (1 - Tbar) * nbar
    J <- (MT/MDES)^2 * numr/denom
  } else if (design == 'd3.3_m3rc2rc')
  {
    numr <- nbar * ICC.2 * (1 - R2.2) + (1 - ICC.2 - ICC.3) * (1 - R2.1)
    denom <- nbar * ( Tbar * (1 - Tbar) * K * (MDES/MT)^2 - ICC.3 * (1 - R2.3) )
    J <- numr/denom
  } else if (design == 'd3.2_m3ff2rc')
  {
    numr <- nbar * ICC.2 * (1 - R2.2) + (1 - ICC.2 - ICC.3) * (1 - R2.1)
    denom <- nbar * Tbar * (1 - Tbar) * K * (MDES/MT)^2
    J <- numr/denom
  } else if (design == 'd3.2_m3rr2rc')
  {
    numr <- nbar * ICC.2 * (1 - R2.2) + (1 - ICC.2 - ICC.3) * (1 - R2.1)
    denom <- nbar * Tbar * (1 - Tbar) * ( K * (MDES/MT)^2 - ICC.3 * omega.3 )
    J <- numr/denom
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
#' @export
calc.K <- function(design, MT, MDES, J, nbar, Tbar,
                   R2.1, R2.2, R2.3,
                   ICC.2, ICC.3,
                   omega.2, omega.3) {

  K <- NA
  if(design == 'd3.1_m3rr2rr')
  {
    K <- (MT/MDES)^2 * ( (ICC.3 * omega.3) +
                           (ICC.2 * omega.2) / J +
                           ((1 - ICC.2 - ICC.3) * (1 - R2.1))/(Tbar * (1 - Tbar) * J * nbar) )
  } else if (design == 'd3.3_m3rc2rc')
  {
    K <- (MT/MDES)^2 * ( (ICC.3 * (1 - R2.3)) / (Tbar * (1 - Tbar)) +
                           (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
                           ((1 - ICC.2 - ICC.3)*(1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
  } else if (design == 'd3.2_m3ff2rc')
  {
    K <- (MT/MDES)^2 * ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
                           ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
  } else if (design == 'd3.2_m3rr2rc')
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



calc_MT = function( df, alpha, two.tailed, target.power ) {
  # t statistics
  T1 <- ifelse(two.tailed == TRUE, abs(qt(alpha/2, df)), abs(qt(alpha, df)))
  T2 <- abs(qt(target.power, df))

  # number of SEs we need for MDES
  MT <- ifelse(target.power >= 0.5, T1 + T2, T1 - T2)

  return(MT)
}



#' Calculating Needed Sample Size for Raw (Unadjusted) Power
#'
#' This is a Helper function for getting a needed Sample Size when no
#' adjustments has been made to the test statistics.
#'
#' It requires iteration because we do not know the degrees of freedom, and so
#' we guess and then calculate sample size, and then recalculate df based on
#' sample size, until we converge.
#'
#' It is possible that the returned sample size will be the minimum sample size
#' required to have at least 1 degree of freedom (even if this provides higher
#' than target level power).
#'
#' @inheritParams pump_power
#'
#' @param design a single RCT design (see list/naming convention)
#' @param typesample type of sample size to calculate: J, K, or nbar
#' @param target.power target power to arrive at
#' @param two.tailed whether to calculate two-tailed or one-tailed power
#' @param max.steps how many steps allowed before terminating
#' @param warn.small Warn if degrees of freedom issues are causing inability to
#'   achieve target power for sample size.
#'
#' @return Requisit sample size (as integer) and associated degreess of freedom.
#'
#' @export
pump_sample_raw <- function(
  design, MTP, typesample,
  MDES,
  nbar = NULL, J = NULL, K = NULL,
  target.power,
  Tbar, alpha, two.tailed = TRUE,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2 = NULL, ICC.3 = NULL,
  omega.2 = NULL, omega.3 = NULL, max.steps = 100,
  warn.small = FALSE
)
{
  if ( typesample == "nbar" ) {
    stopifnot( is.null( nbar ) )
    nbar <- Inf
  } else if ( typesample == "J" ) {
    stopifnot( is.null( J ) )
    J <- Inf
  } else if ( typesample == "K" ) {
    stopifnot( is.null( K ) )
    K <- Inf
  }

  initial_df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
  stopifnot( initial_df > 0 )

  i <- 0
  conv <- FALSE

  # Get initial size (will be low)
  MT <- calc_MT(df = initial_df, alpha = alpha, two.tailed = two.tailed, target.power = target.power)
  if (typesample == "J") {
    J <- calc.J( design, MT = MT, MDES = MDES[1], K = K, nbar = nbar, Tbar = Tbar,
                 R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
                 ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
                 omega.2 = omega.2, omega.3 = omega.3 )
    J <- round(J)
  } else if (typesample == "K") {
    K <- calc.K(
      design, MT = MT, MDES = MDES[1], J = J, nbar = nbar, Tbar = Tbar,
      R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
      ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
      omega.2 = omega.2, omega.3 = omega.3
    )
    K <- round(K)
  } else if (typesample == "nbar") {
    nbar <- calc.nbar(
      design, MT = MT, MDES = MDES[1], J = J, K = K, Tbar = Tbar,
      R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
      ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
      omega.2 = omega.2, omega.3 = omega.3
    )
  }

  df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3, validate = FALSE)

  if( df < 1 ) {
    while( df < 1 ) {
      if ( typesample=="nbar" ) {
        nbar <- nbar + 1
        min_samp_size <- nbar
      } else if ( typesample == "J" ) {
        J <- J + 1
        min_samp_size <- J
      } else if ( typesample == "K" ) {
        K <- K + 1
        min_samp_size <- K
      }
      df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3, validate = FALSE)
    }
    if ( warn.small ) {
      warning(
        'Nonnegative df requirement driving minimum sample size.
        Current sample size will give overpowered study.'
      )
    }
  }


  # Up sample size until we hit our sweet spot.
  while (i <= max.steps & conv == FALSE) {
    df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
    MT <- calc_MT(df = df, alpha = alpha, two.tailed = two.tailed, target.power = target.power)

    if (typesample == "J") {
      J1 <- calc.J( design, MT = MT, MDES = MDES[1], K = K, nbar = nbar, Tbar = Tbar,
                    R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
                    ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
                    omega.2 = omega.2, omega.3 = omega.3 )
      J1 <- round( J1 )

      #cat( "J=", J, "\tdf=", df, "\tJ1=", J1, "\n" )

      if ( is.na(J1) || (J1 <= J) ) {
        conv <- TRUE
      } else {
        J <- J + 1
      }

    } else if (typesample == "K") {
      K1 <- calc.K(
        design, MT = MT, MDES = MDES[1], J = J, nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
        ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      K1 <- round( K1 )
      if ( is.na(K1) || (K1 <= K) ) {
        conv <- TRUE
      } else {
        K <- K + 1
      }
    } else if (typesample == "nbar") {
      nbar1 <- calc.nbar(
        design, MT = MT, MDES = MDES[1], J = J, K = K, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
        ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      #cat( "nbar=", nbar, "\tdf=", df, "\tnbar1=", nbar1, "\n" )

      if (is.na( nbar1 ) || (nbar1 <= nbar) ) {
        conv <- TRUE
      } else {
        nbar <- nbar + 1
      }
    }

    i <- i + 1
  }

  if ( i >= max.steps ) {
    error( "Hit maximum iterations in pump_sample_raw()" )
  }

  if (typesample == "J") {
    if(!is.na(J) & J <= 0){ J <- NA }
    return( list( J=J, df=df ) )
  } else if (typesample == "K") {
    if(!is.na(K) & K <= 0){ K <- NA }
    return( list( K = K, df=df ) )
  } else if (typesample == "nbar") {
    if(!is.na(nbar) & nbar <= 0){ nbar <- NA }
    return( list( nbar=nbar, df=df ) )
  }
}






#' Calculate sample size
#'
#' Given other design parameters, do a search to find the needed sample size to
#' achieve given target level of power to detect specified effect size.
#'
#' @inheritParams pump_mdes
#'
#' @param typesample type of sample size to calculate: "nbar", "J", or "K".
#' @param MDES scalar, or vector of length M; the MDES values for each outcome.
#'
#' @return sample size results
#' @export

pump_sample <- function(
  design, MTP = NULL, typesample,
  MDES, M, numZero = NULL,
  nbar = NULL, J = NULL, K = NULL,
  target.power, power.definition,
  alpha, two.tailed = TRUE,
  Tbar,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  rho = NULL, rho.matrix = NULL,
  omega.2 = 0, omega.3 = 0,
  tnum = 10000, B = 1000,
  max.steps = 20, max.tnum = 2000, start.tnum = 1000, final.tnum = 4*max.tnum,
  cl = NULL, updateProgress = NULL,
  max_sample_size_nbar = 10000,
  max_sample_size_JK = 1000,
  tol = 0.01, give.optimizer.warnings = FALSE,
  just.result.table = TRUE,
  use.logit = FALSE,
  verbose = FALSE )  {

  if ( verbose ) {
    scat( "pump_mdes with %d max iterations per search, starting at %d iterations with final %d iterations.\n\tMax steps %d\n\t%d perms for WY if used\n",
          max.tnum, start.tnum, final.tnum, max.steps, B )
  }


  # Give prelim values for the validation of parameters process.
  if ( typesample == "nbar" ) {
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
    MTP = MTP, MDES = MDES, M = M, J = J, K = K, numZero = numZero,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho, rho.matrix = rho.matrix, B = B
  )
  ##
  params.list <- validate_inputs(design, params.list, single.MDES = TRUE)
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

  # power definition type
  pdef <- parse_power_definition( power.definition, M )

  # validate MTP
  if(MTP == 'None' & !pdef$indiv )
  {
    stop('For minimum or complete power, you must provide a MTP.')
  }

  # Delete parameter we are actually going to search over.
  if ( typesample == "nbar" ) {
    nbar <- NULL
  } else if ( typesample == "J" ) {
    J <- NULL
  } else if ( typesample == "K" ) {
    K <- NULL
  }

  output.colnames <- c("MTP", "Sample type", "Sample size",
                       paste(power.definition, "power") )

  # check if zero power
  if(round(target.power, 2) == 0)
  {
    message('Target power of 0 requested')
    ss.results <- data.frame(MTP, typesample, 0, 0)
    colnames(ss.results) <- output.colnames
    return( make.pumpresult( ss.results, type="sample", params.list=params.list,
                             tries = NULL,
                             just.result.table = just.result.table ) )
  }

  # Checks on what we are estimating, sample size
  if ( verbose ) {
    message(paste("Estimating sample size of type", typesample, "for",
                  MTP, "for target",
                  power.definition, "power of", round(target.power, 4)))
  }

  # Progress Message for the Type of Sample we are estimating, the type of power
  # and the targeted power value
  if (is.function(updateProgress)) {
    msg <- (paste("Estimating", whichSS, "for target", power.definition,
                  "power of",round(power,4)))
    updateProgress(message = msg)
  }

  # adjust bounds to capture needed range
  # for minimum or complete power, expand bounds
  # note: complete power is a special case of minimum power
  if ( !pdef$min ) {
    # Compute needed sample size for raw and BF SS for INDIVIDUAL POWER. We are
    # estimating (potential) bounds
    ss.low <- pump_sample_raw(
      design = design, MTP = MTP, typesample = typesample,
      MDES = MDES, J = J, K = K,
      target.power = target.power,
      nbar = nbar, Tbar = Tbar,
      alpha = alpha,
      two.tailed = two.tailed,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      omega.2 = omega.2, omega.3 = omega.3,
      warn.small = FALSE)

    # Identify sample size for Bonferroni
    ss.high <- pump_sample_raw(
      design = design, MTP = MTP, typesample = typesample,
      MDES = MDES, J = J, K = K,
      target.power = target.power,
      nbar = nbar, Tbar = Tbar,
      alpha = alpha / M, # adjust alpha for BF
      two.tailed = two.tailed,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      omega.2 = omega.2, omega.3 = omega.3,
      warn.small = FALSE)

  } else {
    # lower bound needs to be lower for min type power
    need_pow <- 1 - (1 - target.power)^(1/M)
    ss.low <- pump_sample_raw(
      design = design, MTP = MTP, typesample = typesample,
      MDES = MDES, J = J, K = K,
      target.power = need_pow,
      nbar = nbar, Tbar = Tbar,
      alpha = alpha,
      two.tailed = two.tailed,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      omega.2 = omega.2, omega.3 = omega.3,
      warn.small = FALSE)

    # higher bound needs to be higher for min type power (including comlpete)
    need_pow <- (target.power^(1/M))
    ss.high <- pump_sample_raw(
      design = design, MTP = MTP, typesample = typesample,
      MDES = MDES, J = J, K = K,
      target.power = need_pow,
      nbar = nbar, Tbar = Tbar,
      alpha = alpha / M, # adjust alpha for BF
      two.tailed = two.tailed,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      omega.2 = omega.2, omega.3 = omega.3,
      warn.small = FALSE)
  }

  ss.low.df <- ss.low$df
  ss.low <- ss.low[[1]]

  ss.high.df <- ss.high$sf
  ss.high <- ss.high[[1]]

  # Done if Bonferroni is what we are looking for
  if (MTP == "Bonferroni" & pdef$indiv ) {
    ss.results <- data.frame(MTP, typesample, ss.high, target.power)
    colnames(ss.results) <- output.colnames
    return( make.pumpresult( ss.results, tries = NULL,
                             type="sample", params.list=params.list,
                             just.result.table = just.result.table ) )
  }

  # If we can't make it work with raw, then we can't make it work.
  if ( is.na( ss.low ) ) {
    # ss.results <- data.frame(MTP, typesample, NA, target.power)
    # colnames(ss.results) <- output.colnames
    # return(ss.results)
    ss.low <- 1
  }

  if ( is.na( ss.high ) ) {
    if( typesample == 'nbar')
    {
      ss.high <- max_sample_size_nbar
    } else
    {
      ss.high <- max_sample_size_JK
    }
    warning( "Using default max sample size for one end of initial bounds of search, so estimation may take more time." )

    # Why is this here?   It shouldn't be?
    # start.tnum <- 2000
  }

  ss.low <- round( ss.low )
  ss.high <- round( ss.high )

  # # sometimes we already know the answer!
  # if(ss.low == ss.high)
  # {
  #   ss.results <- data.frame(MTP, typesample, 1, target.power)
  #   colnames(ss.results) <- output.colnames
  #   return(ss.results)
  # }

  # search in the grid from min to max.
  test.pts <- optimize_power(
    design = design, search.type = typesample,
    MTP, target.power, power.definition, tol,
    start.tnum = start.tnum, start.low = ss.low, start.high = ss.high,
    MDES = MDES,
    J = J, K = K, nbar = nbar,
    M = M, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
    rho = rho, omega.2 = omega.2, omega.3 = omega.3,
    B = B, cl = cl,
    max.steps = max.steps, max.tnum = max.tnum,
    final.tnum = final.tnum,
    use.logit = use.logit,
    give.warnings = give.optimizer.warnings,
    verbose = verbose
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

  return( make.pumpresult( ss.results, type = "sample", params.list = params.list,
                           just.result.table = just.result.table,
                           tries = test.pts ) )
}


