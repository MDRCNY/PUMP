
#' List all the supported designs of the `pum` package.
#'
#' List all supported designs, with brief descriptions.
#'
#' @export
supported_designs <- function() {
  cat( "Supported designs:\n" )
  cat( "blocked_i1_2c - Individual rand at level 1, constant impact model\n")
  cat( "blocked_i1_2f - Individual rand at level 1, fixed effects, constant impact model\n")
  cat( "blocked_i1_2r - Individual rand at level 1, random effect for impact (RIRC)\n")
  cat( "blocked_i1_3r - \n")
  cat( "simple_c2_2r - \n")
  cat( "simple_c3_3r - \n")
  cat( "blocked_c2_3f - Randomization at level 2, fixed effects for level 3\n")
  cat( "blocked_c2_3r - Randomization at level 2, random effects\n")

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
    Q.m <- sqrt( (1 - R2.1) /(Tbar * (1-Tbar) * J * nbar) )
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
                 ( ((1 - ICC.2) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * K * nbar) ) )
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
#' @param design a single RCT design (see list/naming convention)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not including block
#'   dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#'
#' @return the degree of freedom

calc.df <- function(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3) {

  if(design == 'blocked_i1_2c')
  {
    df <- J * nbar - numCovar.1 - J - 1
  } else if (design == 'blocked_i1_2f')
  {
    df <- J * nbar - numCovar.1 - 2 * J
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

#' Calculates J, the number of schools
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MT multiplier
#' @param MDES scalar, or vector of length M; the MDES values for each outcome
#' @param J scalar; the number of schools
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects variability
#'
#' @return J, the number of schools

calc.J <- function(design, MT, MDES, nbar, Tbar, R2.1, R2.2, ICC.2, omega.2) {

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
calc.K <- function(design, MT, MDES, J, nbar, Tbar, R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3) {

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


#' Calculate power using PUMP method
#'
#' This functions calculates power for all definitions of power (individual, d-minimal, complete) for all the different MTPs
#' (Bonferroni, Holms, Bejamini-Hocheberg, Westfall-Young Single Step, Westfall-Young Step Down).
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported options: Bonferroni, BH,
#'   Holm, WY-SS, WY-SD
#' @param params.list
#'
#' @return params.list
#' @export
#'
validate_inputs <- function(
  design, MTP, params.list
)
{

  if ( length( MTP ) > 1 ) {
    stop( 'Please provide only a single MTP procedure.' )
  }

  if(length(params.list$MDES) < params.list$M)
  {
    if ( length(params.list$MDES) == 1 ) {
      params.list$MDES <- rep( params.list$MDES, params.list$M )
      warning('Assuming same MDES for all outcomes.  Specify full vector to remove this message.')
    } else {
      stop(paste('Please provide a vector of MDES values of length 1 or M. Current vector:', MDES, 'M =', M))
    }
  }

  if(params.list$K <= 0 | params.list$J <= 0 | params.list$nbar <= 0)
  {
    stop('Please provide positive values of J, K, nbar')
  }

  if(params.list$numCovar.1 < 0 | params.list$numCovar.2 < 0  | params.list$numCovar.3 < 0 )
  {
    stop('Please provide non-negative values of num.Covar')
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

  if(params.list$omega.2 < 0 | params.list$omega.3 < 0)
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
    if(!is.null(params.list$K) | !is.null(params.list$numCovar.3) | !is.null(params.list$R2.3))
    {
      warning('The following parameters are not valid for two-level designs, and will be ignored: K, numCovar.3, R2.3, ICC.3, omega.3')
      params.list$K <- NULL
      params.list$numCovar.3 <- NULL
      params.list$R2.3 <- NULL
      params.list$ICC.3 <- NULL
      params.list$omega.3 <- NULL
    }
  }

  # three level models
  # two level models
  if(design %in% c('blocked_i1_3r', 'simple_c3_3r', 'blocked_c2_3f', 'blocked_c2_3r'))
  {
    if(is.null(params.list$K) | is.null(params.list$numCovar.3) | !is.null(params.list$R2.3))
    {
      stop('The following parameters are required for three-level designs: K, numCovar.3, R2.3')
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
      warning('Omega is assumed to be 0 for constant treatment effects models')
      params.list$omega.2 <- 0
    }
    if(!is.null(params.list$omega.3) && params.list$omega.3 > 0)
    {
      warning('Omega is assumed to be 0 for constant treatment effects models')
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
      warning('Omega2 is assumed to be 0 for blocked_c2_3f model')
      params.list$omega.2 <- 0
    }
    if(params.list$R2.3 > 0)
    {
      warning('R2.3 is assumed to be 0 for blocked_c2_3f model')
      params.list$R2.3 <- 0
    }
    if(params.list$ICC.3 > 0)
    {
      warning('ICC.3 is assumed to be 0 for blocked_c2_3f model')
      params.list$ICC.3 <- 0
    }
  }

  return(params.list)

}

#' Calculate power using PUMP method
#'
#' This functions calculates power for all definitions of power (individual, d-minimal, complete) for all the different MTPs
#' (Bonferroni, Holms, Bejamini-Hocheberg, Westfall-Young Single Step, Westfall-Young Step Down).
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported options: Bonferroni, BH,
#'   Holm, WY-SS, WY-SD
#' @param MDES scalar, or vector of length M; the MDES values for each outcome.
#' @param M scalar; the number of hypothesis tests (outcomes)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param alpha scalar; the family wise error rate (FWER)
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not including block
#'   dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects
#'   variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects
#'   variability
#' @param rho scalar; correlation between outcomes
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
  design, MTP, MDES, M,
  J, K = NULL, nbar, Tbar,
  alpha = 0.05,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1, R2.2 = NULL, R2.3 = NULL,
  ICC.2, ICC.3 = NULL,
  omega.2, omega.3 = NULL,
  rho,
  tnum = 10000, B = 1000,
  cl = NULL,
  updateProgress = NULL
)
{
  # validate input parameters
  params.list = list(
    MDES = MDES, M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho
  )
  ##
  params.list <- validate_inputs(design, MTP, params.list)
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
  rawt.matrix <- mvtnorm::rmvt(tnum, sigma = sigma, df = t.df) + t.shift.mat
  rawp.matrix <- pt(-abs(rawt.matrix), df = t.df) * 2

  if (is.function(updateProgress) & !is.null(rawp.matrix)) {
    updateProgress(message = "P-values have been generated!")
  }

  grab.pval <- function(...,proc) {return(...$adjp[order(...$index),proc])}

  if (MTP == "Bonferroni"){

    adjp <- apply(rawp.matrix, 1, multtest::mt.rawp2adjp, proc = "Bonferroni", alpha = alpha)
    adjp <- do.call(rbind, lapply(adjp, grab.pval, proc = "Bonferroni"))

  } else if (MTP == "Holm") {

    adjp <- apply(rawp.matrix, 1, multtest::mt.rawp2adjp, proc = "Holm", alpha = alpha)
    adjp <- do.call(rbind, lapply(adjp, grab.pval, proc = "Holm"))

  } else if (MTP == "BH") {

    adjp <- apply(rawp.matrix, 1, multtest::mt.rawp2adjp, proc = c("BH"), alpha = alpha)
    adjp <- do.call(rbind, lapply(adjp, grab.pval, proc = "BH"))

  } else if(MTP == "rawp") {

    adjp <- rawp.matrix

  } else if (MTP == "WY-SS"){

    adjp <- adjp.wyss(rawt.matrix = rawt.matrix, B = B, sigma = sigma, t.df = t.df)

  } else if (MTP == "WY-SD"){

    adjp <- adjp.wysd(rawt.matrix = rawt.matrix, B = B, sigma = sigma, t.df = t.df, cl = cl)

  } else {

    stop(paste("Unknown MTP:", MTP))
  }

  if (is.function(updateProgress) & !is.null(adjp)){
    updateProgress(message = paste("Multiple adjustments done for", MTP))
  }

  adjp.each <- list(rawp.matrix, adjp)

  # get matrix of indicators for whether the adjusted p-value is less than alpha
  reject <- function(x) { as.matrix(1*(x < alpha)) }
  reject.each <- lapply(adjp.each, reject)

  # for true positive outcomes, count number of p-values less than 0.05,
  lt.alpha <- function(x) { apply(as.matrix(x[,MDES > 0]), 1, sum) }
  lt.alpha.each <- lapply(reject.each, lt.alpha)

  # individual power
  # mean of columns of booleans of whether adjusted p-values were less than alpha
  power.ind.fun <- function(x) { apply(x, 2, mean) }
  power.ind.each <- lapply(reject.each, power.ind.fun)
  power.ind.each.mat <- do.call(rbind, power.ind.each)

  # mean individual power
  mean.ind.power <- apply(as.matrix(power.ind.each.mat[,MDES > 0]), 1, mean)

  if (is.function(updateProgress) & !is.null(power.ind.each.mat)) {
    updateProgress(message = "Individual power calculation is done.")
  }

  # calculate minimum and complete power
  power.min.fun <- function(x, M) {
    power.min <- numeric(M)
    for (m in 1:M) {
      power.min[m] <- mean(x >= m)
    }
    return(power.min)
  }
  power.min <- lapply(lt.alpha.each, power.min.fun, M = M)
  power.min.mat <- do.call(rbind, power.min)

  # combine all power for all definitions
  all.power.results <- cbind(power.ind.each.mat, mean.ind.power, power.min.mat)

  if(M == 1)
  {
    colnames(all.power.results) = c(paste0("D", 1:M, "indiv"), "indiv.mean", "min", "complete")
  } else
  {
    colnames(all.power.results) = c(paste0("D", 1:M, "indiv"), "indiv.mean", paste0("min",1:(M-1)), "complete")
  }
  rownames(all.power.results) <- c("rawp", MTP)

  if (is.function(updateProgress) & !is.null(all.power.results)) {
    updateProgress(message = paste0("All definitions of power calculation are done."))
  }

  return(all.power.results)
}



#' Midpoint function
#'
#' Calculating the midpoint between the lower and upper bound by calculating
#' half the distance between the two and adding the lower bound to it. The
#' function is a helper function in determining the MDES that falls within
#' acceptable power range.
#'
#' @param lower lower bound
#' @param upper upper bound
#' @importFrom stats dist
#' @return returns midpoint value
midpoint <- function(lower, upper) {
  return(lower + dist(c(lower, upper))[[1]]/2)
}


#' Optimizes power to help in search for MDES or SS
#'
#' @param design a single RCT design (see list/naming convention)
#' @param search.type options: MDES, J, K (nbar not currently supported)
#' @param MTP a single multiple adjustment procedure of interest. Supported options: Bonferroni, BH,
#'   Holm, WY-SS, WY-SD
#' @param target.power Target power to arrive at
#' @param power.definition must be a valid power type outputted by power function, i.e. D1indiv, min1, etc.
#' @param tol tolerance for target power
#' @param start.tnum number of samples for initial power search (should be low to increase speed)
#' @param start.low lower bound for possible power
#' @param start.hight upper bound for possible power
#' @param MDES scalar, or vector of length M; the MDES values for each outcome.
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param alpha scalar; the family wise error rate (FWER)
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not including block
#'   dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects variability
#' @param rho scalar; correlation between outcomes
#' @param tnum scalar; the number of test statistics (samples)
#' @param B scalar; the number of samples/permutations for Westfall-Young
#' @param cl cluster object to use for parallel processing
#' @param max.steps how many steps allowed before terminating
#' @param max.cum.tnum maximum cumulative number of samples
#' @param final.tnum number of samples for final draw
#'
#' @return power
optimize_power <- function(design, search.type, MTP, target.power, power.definition, tol,
                           start.tnum, start.low, start.high,
                           MDES = NULL, J = NULL, K = NULL, nbar = NULL,
                           M = M, Tbar = Tbar, alpha,
                           numCovar.1, numCovar.2, numCovar.3 = NULL,
                           R2.1, R2.2, R2.3 = NULL, ICC.2, ICC.3 = NULL,
                           omega.2, omega.3 = NULL, rho,
                           B = NULL, cl = NULL,
                           max.steps = 20, max.cum.tnum = 5000, final.tnum = 10000)
{

  # fit initial quadratic curve
  # generate a series of points to try
  test.pts <- data.frame(
    step = 0,
    pt = seq(start.low, start.high, length.out = 5),
    power = NA,
    w = start.tnum,
    MTP = MTP,
    target.power = target.power
  )

  # generate power for all these points
  for(i in 1:nrow(test.pts))
  {
    if(search.type == 'mdes'){ MDES <- rep(test.pts$pt[i], M) }
    pt.power.results <- pump_power(
      design, MTP = MTP,
      MDES = MDES,
      J = ifelse(search.type == 'J', test.pts$pt[i], J),
      K = ifelse(search.type == 'K', test.pts$pt[i], K),
      nbar = ifelse(search.type == 'nbar', test.pts$pt[i], nbar),
      tnum = start.tnum,
      # fixed params
      M = M, Tbar = Tbar, alpha = alpha,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      rho = rho, omega.2 = omega.2, omega.3 = omega.3,
      B = B, cl = cl
    )
    if(!(power.definition %in% colnames(pt.power.results)))
    {
      stop(paste0(
        'Please provide a valid power definition. Provided definition: ', power.definition,
        '. Available options: ', paste(colnames(pt.power.results), collapse = ',')))
    }
    test.pts$power[i] <- pt.power.results[MTP, power.definition]
  }

  current.try <- find_best(test.pts, start.low, start.high, target.power, alternate = midpoint(start.low, start.high))
  current.power <- 0
  current.tnum <- start.tnum
  cum.tnum <- 0
  step <- 0

  # fit quadratic based on initial points
  while( (step < max.steps) & (abs( current.power - target.power ) > tol) )
  {
    step <- step + 1
    current.tnum <- pmin(max.cum.tnum, round(current.tnum * 1.1))
    cum.tnum <- cum.tnum + current.tnum

    if(search.type == 'mdes'){ MDES <- rep(current.try, M) }
    current.power.results <- pump_power(
      design, MTP = MTP,
      MDES = MDES,
      J = ifelse(search.type == 'J', current.try, J),
      K = ifelse(search.type == 'K', current.try, K),
      nbar = ifelse(search.type == 'nbar', current.try, nbar),
      tnum = current.tnum,
      # fixed params
      M = M, Tbar = Tbar, alpha = alpha,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      rho = rho, omega.2 = omega.2, omega.3 = omega.3, B = B, cl = cl
    )
    current.power <- current.power.results[MTP, power.definition]

    if(abs(current.power - target.power) < tol)
    {
      check.power.tnum <- pmin(10 * current.tnum, max.cum.tnum)

      if(search.type == 'mdes'){ MDES <- rep(current.try, M) }
      check.power.results <- pump_power(
        design, MTP = MTP,
        MDES = MDES,
        J = ifelse(search.type == 'J', current.try, J),
        K = ifelse(search.type == 'K', current.try, K),
        nbar = ifelse(search.type == 'nbar', current.try, nbar),
        tnum = check.power.tnum,
        # fixed params
        M = M,  Tbar = Tbar, alpha = alpha,
        numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
        R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
        rho = rho, omega.2 = omega.2, omega.3 = omega.3, B = B, cl = cl
      )
      check.power <- check.power.results[MTP, power.definition]

      # cum.tnum <- cum.tnum + check.power.tnum
      # TODO: replace with weighted average?
      current.power <- check.power
      # If still good, go to our final check to see if we are winners!
      # TODO: && (test_pow_R < MAX_ITER)
      if(abs(current.power - target.power) < tol)
      {
        if(search.type == 'mdes'){ MDES <- rep(current.try, M) }
        final.power.results <- pump_power(
          design, MTP = MTP,
          MDES = MDES,
          J = ifelse(search.type == 'J', current.try, J),
          K = ifelse(search.type == 'K', current.try, K),
          nbar = ifelse(search.type == 'nbar', current.try, nbar),
          tnum = final.tnum,
          # fixed params
          M = M, Tbar = Tbar, alpha = alpha,
          numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
          R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
          rho = rho, omega.2 = omega.2, omega.3 = omega.3, B = B, cl = cl
        )
      }
    } # end if within tolerance

    iter.results <- data.frame(
      step = step, pt = current.try, power = current.power, w = current.tnum,
      MTP = MTP, target.power = target.power
    )
    test.pts <- bind_rows(test.pts, iter.results)

    if(current.power < target.power) {
      current.try <- find_best(test.pts, start.low, start.high, target.power, alternate = current.try + 0.10 * (start.high - current.try))
    } else {
      current.try <- find_best(test.pts, start.low, start.high, target.power, alternate = current.try - 0.10 * (current.try - start.low) )
    }
  }

  if( (cum.tnum == max.cum.tnum | step == max.steps) & abs(current.power - target.power) > tol) {
    message("Reached maximum iterations without converging on MDES estimate within tolerance.")
    test.pts <- rbind(test.pts, c(step, NA, NA, NA, MTP, target.power))
  }

  return(test.pts)
}

#' Extract roots from quadratic curve based on given evaluated points
#'
#' @param test.pts power evaluated at different points
#' @param start.low lower bound
#' @param start.high upper bound
#' @param target.power goal power
#' @param alternate alternate point to return if quadratic fit fails
#'
#' @return root of quadratic curve

find_best <- function(test.pts, start.low, start.high, target.power, alternate = NA)
{
  # fit quadratic curve
  quad.mod <- lm( power ~ 1 + pt + I(pt^2), data = test.pts)
  # extract point where it crosses target power
  cc <- rev( coef( quad.mod ) )
  # Using x = [ b pm sqrt( b^2 - 4a(c-y) ) ] / [2a]
  # first check if root exists
  rt.check <- cc[2]^2 - 4 * cc[1] * (cc[3] - target.power)

  if ( rt.check > 0 ) {
    # point to try
    try.pt <- ( -cc[2] + c(-1,1) * sqrt(rt.check) ) / (2 * cc[1] )
    hits <- (start.low <= try.pt) & (try.pt <= start.high)
    if ( sum( hits ) == 1 ) {
      try.pt <- try.pt[hits]
    } else {
      # error
      cat( "Root concerns\n" )
      try.pt <- alternate
    }
  } else {
    cat( "No roots\n" )
    # error
    try.pt <- alternate
  }
  return(try.pt)
}

#' MDES (minimum detectable effect size) function
#'
#' The minimum detectable effect size function calculates the most feasible minimum detectable effect size
#' for a given MTP, power and power definition. The goal is to find the MDES value that satisfies the tolerance
#' set in the parameter in the power value.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported options: Bonferroni, BH,
#'   Holm, WY-SS, WY-SD
#' @param target.power Target power to arrive at
#' @param power.definition must be a valid power type outputted by power function, i.e. D1indiv, min1, etc.
#' @param tol tolerance for target power
#' @param M scalar; the number of hypothesis tests (outcomes)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param alpha scalar; the family wise error rate (FWER)
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not including block
#'   dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects
#'   variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects
#'   variability
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
  nbar, Tbar, alpha, numCovar.1 = 0, numCovar.2 = 0,
  numCovar.3 = 0, R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2, ICC.3 = NULL,
  rho, omega.2, omega.3 = NULL,
  tnum = 10000, B = 1000,
  max.steps = 20, max.cum.tnum = 5000, start.tnum = 200, final.tnum = 10000,
  cl = NULL, updateProgress = NULL
)
{

  # validate input parameters
  params.list = list(
    MDES = MDES, M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho
  )
  ##
  params.list <- validate_inputs(design, MTP, params.list)
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
  if(round(target.power, 2) == 0)
  {
    message('Target power of 0 requested')
    test.pts <- NULL
    mdes.results <- data.frame(MTP, 0, 0)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return(list(mdes.results = mdes.results, test.pts = test.pts))
  }
  # check if max power, then return infinite MDES
  if(round(target.power, 2) == 1)
  {
    message('Target power of 1 requested')
    test.pts <- NULL
    mdes.results <- data.frame(MTP, Inf, 1)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return(list(mdes.results = mdes.results, test.pts = test.pts))
  }

  sigma <- matrix(rho, M, M)
  diag(sigma) <- 1

  message(paste("Estimating MDES for", MTP, "for target", power.definition, "power of", round(target.power, 4)))

  if (MTP == "WY-SD" && B < 1000){
    warning(paste("For the step-down Westfall-Young procedure, it is recommended that sample (B) be at least 1000. Current B:", B))
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
  crit.beta <- ifelse(target.power > 0.5, qt(target.power, df = t.df), qt(1 - target.power, df = t.df))
  mdes.raw  <- ifelse(target.power > 0.5, Q.m * (crit.alpha + crit.beta), Q.m * (crit.alpha - crit.beta))
  mdes.bf   <- ifelse(target.power > 0.5, Q.m * (crit.alphaxM + crit.beta), Q.m * (crit.alphaxM - crit.beta))

  ### raw or bonferroni ###
  if (MTP == "rawp"){
    mdes.results <- data.frame(MTP, mdes.raw, target.power)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return (list(mdes.results = mdes.results, tries = NULL))
  } else if (MTP == "Bonferroni"){
    mdes.results <- data.frame(MTP, mdes.bf, target.power)
    colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))
    return(list(mdes.results = mdes.results, tries = NULL))
  }

  # MDES will be between MDES.raw and MDES.BF
  mdes.low <- mdes.raw
  mdes.high <- mdes.bf

  test.pts <- optimize_power(design, search.type = 'mdes', MTP, target.power, power.definition, tol,
                             start.tnum, start.low = mdes.low, start.high = mdes.high,
                             MDES = NULL, J = J, K = K, nbar = nbar,
                             M = M, Tbar = Tbar, alpha = alpha,
                             numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
                             R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
                             rho = rho, omega.2 = omega.2, omega.3 = omega.3,
                             B = B, cl = cl,
                             max.steps = max.steps, max.cum.tnum = max.cum.tnum,
                             final.tnum = final.tnum)
  mdes.results <- data.frame(MTP, test.pts$pt[nrow(test.pts)], test.pts$power[nrow(test.pts)])
  colnames(mdes.results) <- c("MTP", "Adjusted MDES", paste(power.definition, "power"))

  return(list(mdes.results = mdes.results, test.pts = test.pts))

} # MDES

#' Calculating Sample for Raw (Unadjusted)
#'
#' This is a Helper function for getting Sample Size when no adjustments has been made to the test statistics.
#' The function starts with PowerUp package function mrss.bira2cl but that function seems to have a bug -
#' Only works if we pass in numeric values and not if we pass in objects that hold those values.
#' Additionally, mrss.bira2cl only computes J, not nbar.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported options: Bonferroni, BH,
#'   Holm, WY-SS, WY-SD
#' @param typesample type of sample size to calculate: J, K, or nbar
#' @param MDES scalar, or vector of length M; the MDES values for each outcome
#' @param target.power target power to arrive at
#' @param tol tolerance
#' @param M scalar; the number of hypothesis tests (outcomes)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param J0 scalar; starting point for J
#' @param K0 scalar; starting point for K
#' @param nbar0 scalar; starting point for nbar
#' @param alpha scalar; the family wise error rate (FWER)
#' @param two.tailed whether to calculate two-tailed or one-tailed power
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not including block
#'   dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects
#'   variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects
#'   variability
#' @param rho scalar; correlation between outcomes
#' @param max.steps how many steps allowed before terminating
#' @return raw sample returns
#' @export

pump_sample_raw <- function(
  design, MTP, typesample,
  MDES, M, J = NULL, K = 1,
  J0 = 10, K0 = 4, nbar0 = 10,
  target.power,
  nbar = NULL, Tbar, alpha, two.tailed = TRUE,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2 = NULL, ICC.3 = NULL,
  rho = NULL, omega.2 = NULL, omega.3 = NULL,
  tol = 0.1, max.steps = 100
)
{

  i <- 0
  # convergence
  conv <- FALSE

  while (i <= max.steps & conv == FALSE) {
    # checking which type of sample we are estimating
    if (typesample == "J"){
      df <- calc.df(design, J0, K, nbar, numCovar.1, numCovar.2, numCovar.3)
    } else if (typesample == "K"){
      df <- calc.df(design, J, K0, nbar, numCovar.1, numCovar.2, numCovar.3)
    } else if (typesample == "nbar") {
      df <- calc.df(design, J, K, nbar0, numCovar.1, numCovar.2, numCovar.3)
    }

    # t statistics
    T1 <- ifelse(two.tailed == TRUE, abs(qt(alpha/2, df)), abs(qt(alpha, df)))
    T2 <- abs(qt(target.power, df))
    # multiplier
    MT <- ifelse(target.power >= 0.5, T1 + T2, T1 - T2)

    if (typesample == "J") {
      J1 <- calc.J(
        design, MT = MT, MDES = MDES[1], nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], ICC.2 = ICC.2[1], omega.2 = omega.2
      )
      if (abs(J1 - J0) < tol) {
        conv <- TRUE
      }
      J0 <- (J1 + J0)/2
    } else if (typesample == "K") {
      K1 <- calc.K(
        design, MT = MT, MDES = MDES[1], J = J, nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1], ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      if (abs(K1 - K0) < tol) {
        conv <- TRUE
      }
      K0 <- (K1 + K0)/2
    } else if (typesample == "nbar") {
      nbar1 <- calc.nbar(
        design, MT = MT, MDES = MDES[1], J = J, K = K, Tbar = Tbar,
        R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
        ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
        omega.2 = omega.2, omega.3 = omega.3
      )
      if (abs(nbar1 - nbar0) < tol) {
        conv <- TRUE
      }
      nbar0 <- (nbar1 + nbar0)/2
    }
    i <- i + 1
  }

  if(df < 0 | is.infinite(df)) {
    message('Problem with starting values resulting in impossible df')
  }

  if (typesample == "J") {
    J <- ifelse(df > 0, ceiling(J1), NA)
    return(J)
  } else if (typesample == "K") {
    K <- ifelse(df > 0, ceiling(K1), NA)
    return(K)
  } else if (typesample == "nbar") {
    nbar <- ifelse(df > 0, ceiling(nbar1), NA)
    return(nbar)
  }
}

#' Calculate sample size
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP a single multiple adjustment procedure of interest. Supported options: Bonferroni, BH,
#'   Holm, WY-SS, WY-SD
#' @param typesample type of sample size to calculate: J, K, or nbar
#' @param MDES scalar, or vector of length M; the MDES values for each outcome.
#' @param target.power target power to arrive at
#' @param power.definition must be a valid power type output by power function, i.e. D1indiv, min1, etc.
#' @param tol tolerance
#' @param M scalar; the number of hypothesis tests (outcomes)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the treatment
#' @param J0 scalar; starting point for J
#' @param K0 scalar; starting point for K
#' @param nbar0 scalar; starting point for nbar
#' @param alpha scalar; the family wise error rate (FWER)
#' @param two.tailed whether to calculate two-tailed or one-tailed power
#' @param numCovar.1 scalar; number of Level 1 (individual) covariates (not including block
#'   dummies)
#' @param numCovar.2 scalar; number of Level 2 (school) covariates
#' @param numCovar.3 scalar; number of Level 3 (district) covariates
#' @param R2.1 scalar, or vector of length M; percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect size variability to random effects
#'   variability
#' @param omega.3 scalar; ratio of district effect size variability to random effects
#'   variability
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
  MDES, M, J, K = 1, nbar, Tbar,
  J0 = 10, K0 = 4, nbar0 = 10,
  target.power, power.definition, tol,
  alpha, two.tailed = TRUE,
  numCovar.1 = 0, numCovar.2 = 0,
  numCovar.3 = 0, R2.1, R2.2 = NULL, R2.3 = NULL, ICC.2, ICC.3 = NULL,
  rho, omega.2, omega.3 = NULL,
  tnum = 10000, B = 1000,
  max.steps = 20, max.cum.tnum = 5000, start.tnum = 200, final.tnum = 10000,
  cl = NULL, updateProgress = NULL
)
{
  # extra validation
  if(length(MDES) > 1 & length(unique(MDES)) > 1)
  {
    stop('Procedure assumes MDES is the same for all outcomes.')
  }

  # validate input parameters
  params.list = list(
    MDES = MDES, M = M, J = J, K = K,
    nbar = nbar, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
    ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
    rho = rho
  )
  ##
  params.list <- validate_inputs(design, MTP, params.list)
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

  # save out target sample size
  if(typesample == 'J'){
    target.ss <- J
  } else if(typesample == 'K')
  {
    target.ss <- K
  } else if(typesample == 'nbar')
  {
    target.ss <- nbar
  }
  output.colnames <- c("MTP", "Sample type", "Sample size", paste(power.definition, "power"), "Target sample size")

  # check if zero power, then return 0 MDES
  if(round(target.power, 2) == 0)
  {
    message('Target power of 0 requested')
    test.pts <- NULL
    ss.results <- data.frame(MTP, typesample, 0, 0, target.ss)
    colnames(ss.results) <- output.colnames
    return(list(ss.results = ss.results, test.pts = test.pts))
  }

  message(paste(
    "Estimating sample size of type", typesample, "for", MTP, "for target", power.definition,
    "power of", round(target.power, 4))
  )

  if(typesample == "J"){
    J <- NULL
    nbar0 <- NULL
  } else if (typesample == "K") {
    K <- NULL
    J0 <- NULL
    nbar0 <- NULL
  } else if (typesample == "nbar") {
    nbar <- NULL
    J0 <- NULL
    K0 <- NULL
  }

  if (is.function(updateProgress)) {
    msg <- (paste("Estimating", whichSS, "for target", power.definition, "power of",round(power,4)))
    updateProgress(message = msg)
  }

  # Compute raw and BF SS for INDIVIDUAL POWER.
  # We are estimating bounds like we estimated MDES bounds.
  # for now assuming only two tailed tests
  ss.raw <- pump_sample_raw(
    design = design, MTP, typesample,
    MDES, M, J, K,
    J0, K0, nbar0,
    target.power,
    nbar, Tbar, alpha, two.tailed,
    numCovar.1, numCovar.2, numCovar.3,
    R2.1, R2.2, R2.3, ICC.2, ICC.3,
    rho, omega.2, omega.3
  )
  ss.BF <- pump_sample_raw(
    design = design, MTP, typesample,
    MDES, M, J, K,
    J0, K0, nbar0,
    target.power,
    # change alpha for BF
    nbar, Tbar, alpha/M, two.tailed,
    numCovar.1, numCovar.2, numCovar.3,
    R2.1, R2.2, R2.3, ICC.2, ICC.3,
    rho, omega.2, omega.3
  )

  if (MTP == "rawp"){
    raw.ss <- data.frame(MTP, power.definition, ss.raw, typesample, target.power, target.ss)
    colnames(raw.ss) <- output.colnames
    return(raw.ss)
  } else if (MTP == "Bonferroni") {
    ss.BF <- data.frame(MTP, power.definition, ss.BF, typesample, target.power, target.ss)
    colnames(ss.BF) <- output.colnames
    return(ss.BF)
  }

  # Like the MDES calculation, the sample size would be between raw and Bonferroni.
  # There is no adjustment and there is very conservative adjustment
  ss.low <- ss.raw
  ss.high <- ss.BF

  # sometimes we already know the answer!
  if(ss.low == ss.high)
  {
    test.pts <- NULL
    ss.results <- data.frame(MTP, typesample, 1, target.power, target.ss)
    colnames(ss.results) <- output.colnames
    return(list(ss.results = ss.results, test.pts = test.pts))
  }

  test.pts <- optimize_power(
    design = design, search.type = typesample,
    MTP, target.power, power.definition, tol,
    start.tnum, start.low = ss.low, start.high = ss.high,
    MDES = MDES,
    J = ifelse(typesample == 'J', NULL, J),
    K = ifelse(typesample == 'K', NULL, K),
    nbar = ifelse(typesample == 'nbar', NULL, nbar),
    M = M, Tbar = Tbar, alpha = alpha,
    numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
    R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
    rho = rho, omega.2 = omega.2, omega.3 = omega.3,
    B = B, cl = cl,
    max.steps = max.steps, max.cum.tnum = max.cum.tnum,
    final.tnum = final.tnum
  )
  ss.results <- data.frame(
    MTP, typesample, ifelse(is.na(test.pts$pt[nrow(test.pts)]), NA, ceiling(test.pts$pt[nrow(test.pts)])),
    test.pts$power[nrow(test.pts)], target.ss
  )
  colnames(ss.results) <- output.colnames

  return(list(ss.results = ss.results, test.pts = test.pts))
}
