
#' List all the supported designs of the `pum` package.
#'
#' List all supported designs, with brief descriptions.
#'
#' @export
supported_designs <- function() {
  design = tibble::tribble( ~ Code, ~ Comment,
                   "blocked_i1_2c", "Individual rand at level 1, constant impact model",
                   "blocked_i1_2f", "Individual rand at level 1, fixed effects, constant impact model",
                   "blocked_i1_2r", "Individual rand at level 1, random effect for impact (RIRC)",
                   "blocked_i1_3r", "",
                   "simple_c2_2r", "",
                   "simple_c3_3r", "",
                   "blocked_c2_3f", "Randomization at level 2, fixed effects for level 3",
                   "blocked_c2_3r", "Randomization at level 2, random effects" )

  adjust = tibble::tribble( ~ Code, ~ Comment,
                            "Bonferroni", "The classic (and conservative) multiple testing correction",
                            "Holm", "Bonferroni improved!",
                            "BH", "Benjamini-Hochberg (False Discovery Rate)",
                            "WY-SS", "Westfall-Young, Single Step",
                            "WY-SD", "Westfall-Young, Step Down" )

  dplyr::bind_rows( Design=design, Adjustment=adjust, .id="Group")
}


#' Computes Q_m, the standard error of the effect size estimate
#'
#' @param design a single RCT design (see list/naming convention)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param R2.1 vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects variability
#'
#' @return Q_m, the standard error of the effect size estimate

calc.Q.m <- function(design, J, K, nbar, Tbar, R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3) {

  if(design %in% c('blocked_i1_2c', 'blocked_i1_2f'))
  {
    Q.m <- sqrt( ( (1 - ICC.2)*(1 - R2.1) ) /(Tbar * (1-Tbar) * J * nbar) )
  } else if (design == 'blocked_i1_2r')
  {
    Q.m <- sqrt( (ICC.2 * omega.2)/J +
                ((1 - ICC.2) * (1 - R2.1)) / (Tbar * (1-Tbar) * J * nbar) )
  } else if (design == 'blocked_i1_3r')
  {
    Q.m <- sqrt( (ICC.3 * omega.3) / K +
                 (ICC.2 * omega.2) / (J * K) +
                 ((1 - ICC.2 - ICC.3) * (1 - R2.1))/(Tbar * (1-Tbar) * J * K * nbar) )
  } else if (design == 'simple_c2_2r')
  {
    Q.m <- sqrt( (ICC.2 * (1 - R2.2)) / (Tbar * (1-Tbar) * J) +
                 (1 - ICC.2)*(1 - R2.1) / (Tbar * (1-Tbar) * J * nbar))
  } else if (design == 'simple_c3_3r')
  {
    Q.m <- sqrt( (ICC.3 * (1 - R2.3)) / (Tbar * (1-Tbar) * K) +
                 (ICC.2 * (1 - R2.2)) / (Tbar * (1-Tbar) * J * K) +
                 ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1-Tbar) * J * K * nbar) )
  } else if (design == 'blocked_c2_3f')
  {
    Q.m <- sqrt( ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K) ) +
                 ( ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * K * nbar) ) )
  } else if (design == 'blocked_c2_3r')
  {
    Q.m <- sqrt( ( (ICC.3 * omega.3) / K ) +
                 ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K) ) +
                 ( ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * K * nbar)))
  } else
  {
    stop(paste('Design not implemented:', design))
  }
  return(Q.m)
}


#' Calculate the degrees of freedom for a particular design
#'
#' Given sample sizes, return the used degrees of freedom (frequently
#' conservative) for the design.
#'
#' @inheritParams pump_power
#'
#' @return Degree of freedom for the design.
#'
#' @export

calc.df <- function(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3) {

  if(design == 'blocked_i1_2c')
  {
    df <- J * (nbar - 1) - numCovar.1 - 1
  } else if (design == 'blocked_i1_2f')
  {
    df <- J * (nbar - 2) - numCovar.1
  } else if (design == 'blocked_i1_2r')
  {
    df <- J - numCovar.1 - 1
  } else if (design == 'blocked_i1_3r')
  {
    df <- K - numCovar.3 - 1
  } else if (design == 'simple_c2_2r')
  {
    df <- J - numCovar.1 - 2
  } else if (design == 'simple_c3_3r')
  {
    df <- K - numCovar.3 - 2
  } else if (design == 'blocked_c2_3f')
  {
    df <- K * (J - 2) - numCovar.2
  }else if (design == 'blocked_c2_3r')
  {
    df <- K - numCovar.3 - 1
  } else
  {
    stop(paste('Design not implemented:', design))
  }
  return(df)
}


#' Validates user inputs
#'
#' This functions takes in a list of user inputs. Depending on the inputs,
#' it produces errors or warnings, and at times modifies inputs if necessary.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest.
#' @param params.list a list of parameters input by a user
#'
#' @return params.list
#'
validate_inputs <- function( design, MTP, params.list,
                             mdes.call = FALSE,
                             single.MDES = FALSE )
{
  if(!(design %in% c('blocked_i1_2c', 'blocked_i1_2f', 'blocked_i1_2r',
                     'blocked_i1_3r', 'simple_c2_2r', 'simple_c3_3r',
                     'blocked_c2_3f', 'blocked_c2_3r')))
  {
    stop('Invalid design.')
  }

  if(length( MTP ) > 1)
  {
    stop( 'Please provide only a single MTP procedure.' )
  }

  if(!(MTP %in% c('rawp', 'Bonferroni', 'BH', 'Holm', 'WY-SS', 'WY-SD')))
  {
    stop('Invalid MTP.')
  }

  if ( mdes.call ) {
    if ( !is.null( params.list$MDES ) ) {
      stop( "You cannot provide MDES to pump_mdes()" )
    }
  } else {
    if(!is.null(params.list$numZero) & length(params.list$MDES) != 1)
    {
      stop('If providing a number of zero outcomes, please provide a single MDES value.')
    }
    if(!is.null(params.list$numZero))
    {
      numNonzero <- params.list$M - params.list$numZero
      params.list$MDES <- c(rep(params.list$MDES, numNonzero), rep(0, params.list$numZero))
      print(paste('Assumed full MDES vector:', params.list$MDES))
    }

    if ( single.MDES ) {
      if ( length(params.list$MDES) != 1 ) {
        stop( "Please provide a single MDES value This function does not support vector MDES inputs." )
      }
    } else if(length(params.list$MDES) < params.list$M)
    {
      if ( length(params.list$MDES) == 1 ) {
        params.list$MDES <- rep( params.list$MDES, params.list$M )
        warning('Assuming same MDES for all outcomes.  Specify full vector to remove this message.')
      } else {
        stop(paste('Please provide a vector of MDES values of length 1 or M. Current vector:',
                   MDES, 'M =', M))
      }
    }

  }

  if( (!is.null( params.list$K ) && params.list$K <= 0) | params.list$J <= 0 | params.list$nbar <= 0)
  {
    stop('Please provide positive values of J, K, and/or nbar')
  }

  if(params.list$numCovar.1 < 0 | params.list$numCovar.2 < 0  |
     ( !is.null( params.list$numCovar.3 ) && params.list$numCovar.3 < 0 ) )
  {
    stop('Please provide non-negative values of your num.Covar parameters')
  }

  if(params.list$Tbar >= 1 | params.list$Tbar <= 0)
  {
    stop('Please provide Tbar as a probability strictly between 0 and 1')
  }

  if(params.list$alpha > 1 | params.list$alpha < 0)
  {
    stop('Please provide alpha as a probability  between 0 and 1')
  }

  if(any(params.list$R2.1 > 1) | any(params.list$R2.1 < 0) |
     any(params.list$R2.2 > 1) | any(params.list$R2.2 < 0) |
     any(params.list$R2.3 > 1) | any(params.list$R2.3 < 0))
  {
    stop('Please provide R2 as a probability between 0 and 1')
  }

  if(params.list$omega.2 < 0 | (!is.null(params.list$omega.3) && params.list$omega.3 < 0))
  {
    stop('Please provide a non-negative value of Omega')
  }

  if(params.list$rho > 1 | params.list$rho < -1)
  {
    stop('Please provide rho as a correlation between -1 and 1')
  }

  # two level models
  if(design %in% c('blocked_i1_2c', 'blocked_i1_2f', 'blocked_i1_2r', 'simple_c2_2r'))
  {
    if( (!is.null(params.list$K) && params.list$K > 1 )|
       ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
       ( !is.null(params.list$R2.3)) && any( params.list$R2.3 > 0 ) |
       ( !is.null(params.list$omega.3) && params.list$omega.3 > 0 ) )
    {
      warning('The following parameters are not valid for two-level designs, and will be ignored: K, numCovar.3, R2.3, ICC.3, omega.3')
      params.list$K <- NULL
      params.list$numCovar.3 <- NULL
      params.list$R2.3 <- NULL
      params.list$ICC.3 <- NULL
      params.list$omega.3 <- NULL
    }
  }

  if(design == 'blocked_c2_3f')
  {
    if( ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
        ( !is.null(params.list$R2.3) && any( params.list$R2.3 > 0 ) ) )
    {
      warning('The following parameters are not valid for fixed effect designs, and will be ignored: numCovar.3, R2.3')
      params.list$numCovar.3 <- NULL
      params.list$R2.3 <- NULL
    }
  }


  # three level models
  if(design %in% c('blocked_i1_3r', 'simple_c3_3r', 'blocked_c2_3f', 'blocked_c2_3r'))
  {
    if(is.null(params.list$K) || params.list$K <= 1 )
    {
      stop('You must specify K, with K > 1 (number of units at level 3) for three-level designs' )
    }
  }

  # three level models, continued.
  if(design %in% c('blocked_i1_3r', 'simple_c3_3r', 'blocked_c2_3r'))
  {
    if( is.null(params.list$numCovar.3) | is.null(params.list$R2.3))
    {
      stop('You must specify both numCovar.3 and R2.3 for three-level designs with random effects')
    }
  }



  # ICC
  if(!is.null(params.list$ICC.2) && !is.null(params.list$ICC.3) && params.list$ICC.2 + params.list$ICC.3 > 1)
  {
    stop('ICC.2 + ICC.3 must be <= 1')
  }

  # constant treatment effects models
  if(design %in% c('blocked_i1_2c', 'simple_c2_2r', 'simple_c3_3r'))
  {
    if(params.list$omega.2 > 0)
    {
      warning('Omega is assumed to be 0 for constant treatment effects models. Ignoring input omega.2 value')
      params.list$omega.2 <- 0
    }
    if(!is.null(params.list$omega.3) && params.list$omega.3 > 0)
    {
      warning('Omega is assumed to be 0 for constant treatment effects models. Ignoring input omega.3 value')
      params.list$omega.3 <- 0
    }
  }

  # specific 3 level models
  if(design %in% c('blocked_i1_3r', 'blocked_c2_3f', 'blocked_c2_3r'))
  {
    if(is.null(params.list$omega.3))
    {
      stop('Omega.3 is required for this design.')
    }
  }
  if(design %in% c('blocked_i1_3r', 'blocked_c2_3r'))
  {
    if(is.null(params.list$ICC.3))
    {
      stop('ICC.3 is required for this design.')
    }
  }

  if(design == 'blocked_c2_3f')
  {
    if(params.list$omega.2 > 0)
    {
      warning('Omega2 is assumed to be 0 for blocked_c2_3f model. Ignoring input omega.2 value')
      params.list$omega.2 <- 0
    }

    # NOTE: I believe we could have a ICC > 0, even if we do not model it.  This
    # would force other variation down. -lwm
    #
    # if(params.list$ICC.3 > 0)
    # {
    #   warning('ICC.3 is assumed to be 0 for blocked_c2_3f model. Ignoring input ICC.3 value')
    #   params.list$ICC.3 <- 0
    # }
  }

  return(params.list)

}


#' Calculates different definitions of power
#'
#' This function takes in a matrix of adjusted p-values and outputs different types of power
#'
#' @param pval.mat matrix of p-values, columns are outcomes
#' @param ind.nonzero vector indicating which outcomes are nonzero
#' @param alpha scalar; the family wise error rate (FWER)
#'
#' @return power results for individual, minimum, complete power
#' @export
get.power.results = function(pval.mat, ind.nonzero, alpha)
{
  M <- ncol(pval.mat)
  num.nonzero <- sum(ind.nonzero)

  # rejected tests
  rejects <- apply(pval.mat, 2, function(x){ 1*(x < alpha) })
  rejects.nonzero <- rejects[,ind.nonzero]

  # individual power
  power.ind <- apply(rejects.nonzero, 2, mean)
  power.ind.mean <- mean(power.ind)

  # minimum power
  power.min <- rep(NA, num.nonzero)
  for(m in 1:num.nonzero)
  {
    min.rejects <- apply(rejects.nonzero, 1, function(x){ sum(x) >= m })
    power.min[m] <- mean(min.rejects)
  }

  # combine all power for all definitions
  all.power.results <- data.frame(matrix(c(power.ind, power.ind.mean, power.min), nrow = 1))
  colnames(all.power.results) = c(paste0("D", 1:num.nonzero, "indiv"), "indiv.mean", paste0("min",1:(num.nonzero-1)), "complete")

  return(all.power.results)
}


#' Calculate power using PUMP method
#'
#' This functions calculates power for all definitions of power (individual,
#' d-minimal, complete) for all the different MTPs (Bonferroni, Holms,
#' Bejamini-Hocheberg, Westfall-Young Single Step, Westfall-Young Step Down).
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP Multiple adjustment procedure of interest. Supported options:
#'   Bonferroni, BH, Holm, WY-SS, WY-SD (passed as strings).  Provide list to
#'   automatically re-run for each procedure on the list.
#' @param MDES scalar, or vector of length M; the MDES values for each outcome.
#' @param numZero if MDES is scalar, number of outcomes assumed to be zero.
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
#'   Level 1 covariates for each outcome. Defaults to 0.
#' @param R2.2 scalar, or vector of length M; percent of variation explained by
#'   Level 2 covariates for each outcome. Defaults to 0.
#' @param R2.3 scalar, or vector of length M; percent of variation explained by
#'   Level 3 covariates for each outcome. Defaults to 0.
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of variance of school-average impacts to
#'   variance of school-level random intercepts.  Default to 0 (no treatment
#'   variation).
#' @param omega.3 scalar; ratio of variance of district-average impacts to
#'   variance of district-level random intercepts. Default to 0 (no treatment
#'   variation).
#' @param rho scalar; assumed correlation between the test statistics of the
#'   tests.
#' @param tnum scalar; the number of test statistics (samples)
#' @param B scalar; the number of samples/permutations for Westfall-Young
#' @param cl cluster object to use for parallel processing
#' @param updateProgress the callback function to update the progress bar (User
#'   does not have to input anything)
#'
#' @importFrom multtest mt.rawp2adjp
#' @return power results for MTP and unadjusted across all definitions of power
#' @export
#'
#'
pump_power <- function(
  design, MTP, MDES, numZero = NULL,
  M, J, K = 1, nbar, Tbar,
  alpha = 0.05,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  omega.2 = 0, omega.3 = 0,
  rho,
  tnum = 10000, B = 3000,
  cl = NULL,
  updateProgress = NULL,
  validate.inputs = TRUE
)
{
  # Call self for each element on MTP list.
  if ( length( MTP ) > 1 ) {
    des = purrr::map( MTP,
                     pum::pump_power, design=design, MDES=MDES, M=M, J=J, K = K, nbar=nbar,
                     Tbar=Tbar,
                     alpha=alpha, numCovar.1 = numCovar.1, numCovar.2 = numCovar.2,
                     numCovar.3 = numCovar.3, R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                     ICC.2 = ICC.2, ICC.3 = ICC.3,
                     rho=rho, omega.2=omega.2, omega.3 = omega.3,
                     tnum = tnum, B = B, cl = cl, updateProgress = updateProgress )
    des = do.call( rbind, des )
    des = des[ -seq( 3, nrow(des), by=2 ), ]
    return( des )
  }

  if(validate.inputs)
  {
    # validate input parameters
    params.list <- list(
      MDES = MDES, numZero = numZero, M = M, J = J, K = K,
      nbar = nbar, Tbar = Tbar, alpha = alpha,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
      ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
      rho = rho
    )

    params.list <- validate_inputs(design, MTP, params.list)

    MDES <- params.list$MDES
    M <- params.list$M; J <- params.list$J; K <- params.list$K
    nbar <- params.list$nbar; Tbar <- params.list$Tbar; alpha <- params.list$alpha
    numCovar.1 <- params.list$numCovar.1; numCovar.2 <- params.list$numCovar.2
    numCovar.3 <- params.list$numCovar.3
    R2.1 <- params.list$R2.1; R2.2 <- params.list$R2.2; R2.3 <- params.list$R2.3
    ICC.2 <- params.list$ICC.2; ICC.3 <- params.list$ICC.3
    omega.2 <- params.list$omega.2; omega.3 <- params.list$omega.3
    rho <- params.list$rho
  }

  # compute test statistics for when null hypothesis is false
  Q.m <- calc.Q.m(
    design = design, J = J, K = K, nbar = nbar, Tbar = Tbar,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3,
    omega.2 = omega.2, omega.3 = omega.3
  )
  t.shift <- MDES/Q.m
  t.df <- calc.df(
    design = design, J = J, K = K,
    nbar = nbar,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3
  )
  t.shift.mat <- t(matrix(rep(t.shift, tnum), M, tnum))

  # correlation between the test statistics
  sigma <- matrix(rho, M, M)
  diag(sigma) <- 1

  # generate t values and p values under alternative hypothesis using multivariate t-distribution
  rawt.mat <- mvtnorm::rmvt(tnum, sigma = sigma, df = t.df) + t.shift.mat
  rawp.mat <- pt(-abs(rawt.mat), df = t.df) * 2

  if (is.function(updateProgress) & !is.null(rawp.mat)) {
    updateProgress(message = "P-values have been generated!")
  }

  grab.pval <- function(...,proc) {return(...$adjp[order(...$index),proc])}

  if (MTP == "Bonferroni"){

    adjp <- apply(rawp.mat, 1, multtest::mt.rawp2adjp, proc = "Bonferroni", alpha = alpha)
    adjp <- do.call(rbind, lapply(adjp, grab.pval, proc = "Bonferroni"))

  } else if (MTP == "Holm") {

    adjp <- apply(rawp.mat, 1, multtest::mt.rawp2adjp, proc = "Holm", alpha = alpha)
    adjp <- do.call(rbind, lapply(adjp, grab.pval, proc = "Holm"))

  } else if (MTP == "BH") {

    adjp <- apply(rawp.mat, 1, multtest::mt.rawp2adjp, proc = c("BH"), alpha = alpha)
    adjp <- do.call(rbind, lapply(adjp, grab.pval, proc = "BH"))

  } else if(MTP == "rawp") {

    adjp <- rawp.mat

  } else if (MTP == "WY-SS"){

    adjp <- adjp.wyss(rawt.mat = rawt.mat, B = B, sigma = sigma, t.df = t.df)

  } else if (MTP == "WY-SD"){

    adjp <- adjp.wysd(rawt.mat = rawt.mat, B = B, sigma = sigma, t.df = t.df, cl = cl)

  } else {

    stop(paste("Unknown MTP:", MTP))
  }

  if (is.function(updateProgress) & !is.null(adjp)){
    updateProgress(message = paste("Multiple adjustments done for", MTP))
  }

  ind.nonzero <- MDES > 0
  power.results.rawp <- get.power.results(rawp.mat, ind.nonzero, alpha)
  power.results.proc <- get.power.results(adjp, ind.nonzero, alpha)
  power.results.all <- data.frame(rbind(power.results.rawp, power.results.proc))
  rownames(power.results.all) <- c('rawp', MTP)

  return(power.results.all)
}





scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}

#' Run pump_power on combination of factors
#'
#' This extenstion of `pump_power()` will take lists of parameter values and run
#' `pump_power()` on all combinations of these values.
#'
#' It can only assume the same MDES value for all outcomes due to this.
#'
#' Each parameter in the parameter list can be a list, not scalar.
#'
#' @inheritParams pump_power
#'
#' @param MDES This is *not* a list of MDES for each outcome, but rather a list
#'   of MDES to explore. Each value will be assumed held constant across all M
#'   outcomes.
#'
#' @export
pump_power_grid <- function( design, MTP, MDES, M, J, K = 1, nbar, Tbar, alpha,
                             numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
                             R2.1, R2.2 = NULL, R2.3 = NULL,
                             ICC.2, ICC.3 = NULL,
                             rho, omega.2 = NULL, omega.3 = NULL, ... ) {

  args = list( M = M, J = J, K = K, nbar = nbar,
               Tbar = Tbar, alpha = alpha,
               numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
               R2.1 = R2.1, R2.2 = R2.2, ICC.2 = ICC.2, ICC.3 = ICC.3,
               rho = rho, omega.2 = omega.2, omega.3 = omega.3 )
  nulls = purrr::map_lgl( args, is.null )
  args = args[ !nulls ]

  grid = do.call( expand.grid, args )
  scat( "Processing %d calls\n", nrow(grid) )
  grid$res = purrr::pmap( grid, pump_power,
                          design = design,
                          MTP = MTP,
                          MDES = MDES, ... )

  grid$res = purrr::map( grid$res, as.data.frame )
  grid$res = purrr::map( grid$res, tibble::rownames_to_column, var ="adjustment" )
  grid = tidyr::unnest( grid, res )

  grid
}
