
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

  if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
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
    denom <- Tbar * (1 - Tbar) * J * ((MDES/MT)^2) - J*ICC.3*(1 - R2.3)  - ICC.2 * (1 - R2.2)
    nbar <- numr / denom
  } else if (design == 'd3.2_m3ff2rc')
  {
    numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J * K * ((MDES/MT)^2) - ICC.2 * (1 - R2.2)
    nbar <- numr / denom
  } else if (design == 'd3.2_m3rr2rc')
  {
    numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
    denom <- Tbar * (1 - Tbar) * J * K * ((MDES/MT)^2) * ICC.3 * omega.3 - ICC.2 * (1 - R2.2)
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

calc.J <- function(design, MT = 2.8, MDES, nbar, Tbar, R2.1, R2.2, ICC.2, omega.2) {

  if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
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
    denom <- Tbar * (1 - Tbar) * K * (MDES/MT)^2 - ICC.3 * (1 - R2.3)
    J <- (1 / nbar) * numr/denom
  } else if (design == 'd3.2_m3ff2rc')
  {
    numr <- nbar * ICC.2 * (1 - R2.2) + (1 - ICC.2 - ICC.3) * (1 - R2.1)
    denom <- nbar * Tbar * (1 - Tbar) * K * (MDES/MT)^2
    J <- numr/denom
  } else if (design == 'd3.2_m3rr2rc')
  {
    numr <- nbar * ICC.2 * (1 - R2.2) + (1 - ICC.2 - ICC.3) * (1 - R2.1)
    denom <- nbar * K * (MDES/MT)^2 - nbar * ICC.3 * omega.3
    J <- 1 / (Tbar * (1 - Tbar)) * numr/denom
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
#'
#' @return Requisit sample size (as integer).
#' @export
pump_sample_raw <- function(
  design, MTP, typesample,
  MDES,
  nbar = NULL, J = NULL, K = NULL,
  target.power,
  Tbar, alpha, two.tailed = TRUE,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2 = NULL, ICC.3 = NULL,
  omega.2 = NULL, omega.3 = NULL, max.steps = 100
)
{
  if ( typesample=="nbar" ) {
    stopifnot( is.null( nbar ) )
    nbar = Inf
  } else if ( typesample == "J" ) {
    stopifnot( is.null( J ) )
    J = Inf
  } else if ( typesample == "K" ) {
    stopifnot( is.null( K ) )
    K = Inf
  }

  initial_df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
  stopifnot( initial_df > 0 )

  i <- 0
  conv <- FALSE

  # Get initial size (will be low)
  MT = calc_MT(df = initial_df, alpha = alpha, two.tailed = two.tailed, target.power = target.power)
  if (typesample == "J") {
    J <- calc.J( design, MT = MT, MDES = MDES[1], nbar = nbar, Tbar = Tbar,
                 R2.1 = R2.1[1], R2.2 = R2.2[1], ICC.2 = ICC.2[1], omega.2 = omega.2 )
    J = round(J)
  } else if (typesample == "K") {
    K <- calc.K(
      design, MT = MT, MDES = MDES[1], J = J, nbar = nbar, Tbar = Tbar,
      R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
      ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
      omega.2 = omega.2, omega.3 = omega.3
    )
    K = round(K)
  } else if (typesample == "nbar") {
    nbar <- calc.nbar(
      design, MT = MT, MDES = MDES[1], J = J, K = K, Tbar = Tbar,
      R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
      ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
      omega.2 = omega.2, omega.3 = omega.3
    )
  }

  df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
  if( df < 1 ) {
    while( df < 1 ) {
      if ( typesample=="nbar" ) {
        nbar = nbar + 1
        min_samp_size = nbar
      } else if ( typesample == "J" ) {
        J = J + 1
        min_samp_size = J
      } else if ( typesample == "K" ) {
        K = K + 1
        min_samp_size = K
      }
      df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
    }
    warning('Nonnegative df requirement driving minimum sample size. Current sample size will give overpowered study.')
  }


  # Up sample size until we hit our sweet spot.
  while (i <= max.steps & conv == FALSE) {
    df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
    MT = calc_MT(df = df, alpha = alpha, two.tailed = two.tailed, target.power = target.power)

    if (typesample == "J") {
      J1 <- calc.J( design, MT = MT, MDES = MDES[1], nbar = nbar, Tbar = Tbar,
                    R2.1 = R2.1[1], R2.2 = R2.2[1], ICC.2 = ICC.2[1], omega.2 = omega.2 )
      J1 = round( J1 )

      #cat( "J=", J, "\tdf=", df, "\tJ1=", J1, "\n" )

      if ( is.na(J1) || (J1 <= J) ) {
        conv <- TRUE
      } else {
        J = J + 1
      }

    } else if (typesample == "K") {
      K1 <- calc.K(
        design, MT = MT, MDES = MDES[1], J = J, nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
        ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      K1 = round( K1 )
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
    warning( "Hit maximum iterations in pump_sample_raw()" )
  }

  if (typesample == "J") {
    return(J)
  } else if (typesample == "K") {
    return(K)
  } else if (typesample == "nbar") {
    return(nbar)
  }
}



# LWM: I dumped binary search because things got confusing.  Perhaps this is the
# better road?  I think I broke it.  :-(

# Old @param tol tolerance of how much our search can change the sample size before
#'   stopping.  (In number of units, so 0.1 is relatively precise.)

pump_sample_raw_old <- function(
  design, MTP, typesample,
  MDES,
  nbar = NULL, J = NULL, K = NULL,
  target.power,
  Tbar, alpha, two.tailed = TRUE,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2 = NULL, ICC.3 = NULL,
  omega.2 = NULL, omega.3 = NULL, max.steps = 100
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

  initial_df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
  stopifnot( initial_df > 0 )

  i <- 0
  # convergence
  conv <- FALSE

  min_samp_size = 0
  bump = 0

  while (i <= max.steps & conv == FALSE) {
    # checking which type of sample we are estimating
    df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
    if( df < 1 ) {
      while( df < 1 ) {
        if ( typesample=="nbar" ) {
          nbar = ceiling( nbar + 1 )
          min_samp_size = nbar
        } else if ( typesample == "J" ) {
          J = ceiling( J + 1 )
          min_samp_size = J
        } else if ( typesample == "K" ) {
          K = ceiling( K + 1 )
          min_samp_size = K
        }
        df <- calc.df(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3)
      }
      warning('Landed on situation with impossible df.  Searching for min sample size to allow estimation with model.')
    }

    # t statistics
    T1 <- ifelse(two.tailed == TRUE, abs(qt(alpha/2, df)), abs(qt(alpha, df)))
    T2 <- abs(qt(target.power, df))

    # number of SEs we need for MDES
    MT <- ifelse(target.power >= 0.5, T1 + T2, T1 - T2)

    #cat( "df = ", as.numeric( df ), "\n" )

    if (typesample == "J") {
      J1 <- calc.J(
        design, MT = MT, MDES = MDES[1], nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], ICC.2 = ICC.2[1], omega.2 = omega.2
      )
      J1 = round( J1 )
      cat( "J=", J, "\tdf=", df, "\tJ1=", J1, "\n" )
      if ( is.na(J1) || (J1 == J) ) {
        conv <- TRUE
      }
      if ( J < min_samp_size ) {
        min_samp_size = min_samp_size + 1
        J = min_samp_size
      }

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
      K <- pmax( round(K1), min_samp_size)
      if ( K == min_samp_size ) {
        min_samp_size = min_samp_size + 1
      }

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
      nbar <- pmax( round(nbar1), min_samp_size )
      if ( nbar == min_samp_size ) {
        min_samp_size = min_samp_size + 1
      }

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
#' @param two.tailed whether to calculate two-tailed or one-tailed power
#' @param rho scalar; correlation between outcomes
#' @param tnum scalar; the number of test statistics (samples)
#' @param B scalar; the number of samples/permutations for Westfall-Young
#' @param max.steps how many steps allowed before terminating
#' @param max.cum.tnum maximum cumulative number of samples
#' @param final.tnum number of samples for final draw
#' @param cl cluster object to use for parallel processing
#' @param tol tolerance of how close power should be to the target power.
#' @param updateProgress the callback function to update the progress bar (User
#'   does not have to input anything)
#' @param max_sample_size If the initial bounds fail to capture a range of
#'   working values (e.g., for non-achievable power on one end of the search and
#'   not the other), default to this sample size in the search.
#'
#' @return sample size results
#' @export

pump_sample <- function(
  design, MTP, typesample,
  MDES, M,
  nbar = NULL, J = NULL, K = NULL,
  target.power, power.definition,
  alpha, two.tailed = TRUE,
  Tbar,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  rho,
  omega.2 = 0, omega.3 = 0,
  tnum = 10000, B = 1000,
  max.steps = 20, max.cum.tnum = 5000, start.tnum = 200, final.tnum = 10000,
  cl = NULL, updateProgress = NULL,
  max_sample_size = 10000,
  tol = 0.01, give.optimizer.warnings = FALSE
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
  params.list <- pum:::validate_inputs(design, MTP, params.list, single.MDES = TRUE)
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

  # Identify sample size for Bonferroni
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

  if ( is.na( ss.high ) ) {
    ss.high <- max_sample_size
    warning( "Using default max sample size for one end of initial bounds of search." )
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
    final.tnum = final.tnum,
    give.warnings = give.optimizer.warnings
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


