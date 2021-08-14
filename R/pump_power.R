
#' List all the supported designs of the `pum` package.
#'
#' List all supported designs, with brief descriptions.
#'
#' @export
supported_designs <- function() {
  design = tibble::tribble( ~ Code, ~ Comment,
                   # 1 level design
                   "d1.1_m2cc", "1 level, level 1 randomization / constant intercepts, constant impacts model",
                   # 2 level designs, randomization at level 1
                   "d2.1_m2fc", "2 lvls, lvl 1 rand / fixed intercepts, constant impacts",
                   "d2.1_m2ff", "2 lvls, lvl 1 rand / fixed intercepts, fixed impacts",
                   "d2.1_m2fr", "2 lvls, lvl 1 rand / fixed intercepts, random impacts",
                   # 3 lvl design, rand at lvl 1
                   "d3.1_m3rr2rr", "3 lvls, lvl 1 rand / lvl 3 random intercepts, random impacts, lvl 2 random intercepts, random impacts",
                   # 2 lvl design, rand at lvl 2
                   "d2.2_m2rc", "2 lvls, lvl 2 rand / random intercepts, constant impacts",
                   # 3 lvl design, rand at lvl 3
                   "d3.3_m3rc2rc", "3 lvls, lvl 3 rand / lvl 3 random intercepts, constant impacts, lvl 2 random intercepts, constant impacts",
                   # 3 lvl design, rand at lvl 2
                   "d3.2_m3ff2rc", "3 lvls, lvl 2 rand / lvl 3 fixed intercepts, fixed impacts, lvl 2 random intercepts, constant impacts",
                   "d3.2_m3rr2rc", "3 lvls, lvl 2 rand / lvl 3 random intercepts, random impacts, lvl 2 random intercepts, constant impacts" )

  design = tidyr::separate( design, Code, into=c("Design","Model"), remove = FALSE, sep="_" )

  adjust = tibble::tribble( ~ Method, ~ Comment,
                            "Bonferroni", "The classic (and conservative) multiple testing correction",
                            "Holm", "Bonferroni improved!",
                            "BH", "Benjamini-Hochberg (False Discovery Rate)",
                            "WY-SS", "Westfall-Young, Single Step",
                            "WY-SD", "Westfall-Young, Step Down" )

  list( Design=design, Adjustment=adjust )
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

  if(design %in% c('d1.1_m2cc'))
  {
    Q.m <- sqrt( ( (1 - R2.1) ) /(Tbar * (1-Tbar) * nbar) )
  } else if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
  {
    Q.m <- sqrt( ( (1 - ICC.2)*(1 - R2.1) ) /(Tbar * (1-Tbar) * J * nbar) )
  } else if (design == 'd2.1_m2fr')
  {
    Q.m <- sqrt( (ICC.2 * omega.2)/J +
                ((1 - ICC.2) * (1 - R2.1)) / (Tbar * (1-Tbar) * J * nbar) )
  } else if (design == 'd3.1_m3rr2rr')
  {
    Q.m <- sqrt( (ICC.3 * omega.3) / K +
                 (ICC.2 * omega.2) / (J * K) +
                 ((1 - ICC.2 - ICC.3) * (1 - R2.1))/(Tbar * (1-Tbar) * J * K * nbar) )
  } else if (design == 'd2.2_m2rc')
  {
    Q.m <- sqrt( (ICC.2 * (1 - R2.2)) / (Tbar * (1-Tbar) * J) +
                 (1 - ICC.2)*(1 - R2.1) / (Tbar * (1-Tbar) * J * nbar))
  } else if (design == 'd3.3_m3rc2rc')
  {
    Q.m <- sqrt( (ICC.3 * (1 - R2.3)) / (Tbar * (1-Tbar) * K) +
                 (ICC.2 * (1 - R2.2)) / (Tbar * (1-Tbar) * J * K) +
                 ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1-Tbar) * J * K * nbar) )
  } else if (design == 'd3.2_m3ff2rc')
  {
    Q.m <- sqrt( ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K) ) +
                 ( ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * K * nbar) ) )
  } else if (design == 'd3.2_m3rr2rc')
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

calc.df <- function(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3, validate = TRUE) {

  if(design == 'd1.1_m2cc')
  {
    df <- J * nbar - numCovar.1 - 1
  } else if(design == 'd2.1_m2fc')
  {
    df <- J * (nbar - 1) - numCovar.1 - 1
  } else if (design == 'd2.1_m2ff')
  {
    df <- J * (nbar - 2) - numCovar.1
  } else if (design == 'd2.1_m2fr')
  {
    df <- J - numCovar.1 - 1
  } else if (design == 'd3.1_m3rr2rr')
  {
    df <- K - numCovar.3 - 1
  } else if (design == 'd2.2_m2rc')
  {
    df <- J - numCovar.1 - 2
  } else if (design == 'd3.3_m3rc2rc')
  {
    df <- K - numCovar.3 - 2
  } else if (design == 'd3.2_m3ff2rc')
  {
    df <- K * (J - 2) - numCovar.2
  }else if (design == 'd3.2_m3rr2rc')
  {
    df <- K - numCovar.3 - 1
  } else
  {
    stop(paste('Design not implemented:', design))
  }

  if(validate & df <= 0)
  {
    stop('Invalid design parameters resulting in nonpositive degrees of freedom')
  }

  return(df)
}


#' Calculates different definitions of power
#'
#' This function takes in a matrix of adjusted p-values and outputs different types of power
#'
#' @param pval.mat matrix of p-values, columns are outcomes
#' @param ind.nonzero vector indicating which outcomes are nonzero
#' @param alpha scalar; the family wise error rate (FWER)
#' @param unadj whether p-values are unadjusted or not
#'
#' @return power results for individual, minimum, complete power
#' @export
get.power.results = function(pval.mat, ind.nonzero, alpha, adj = TRUE)
{
  M <- ncol(pval.mat)
  num.nonzero <- sum(ind.nonzero)

  # rejected tests
  rejects <- apply(pval.mat, 2, function(x){ 1*(x < alpha) })
  rejects.nonzero <- rejects[,ind.nonzero, drop = FALSE]

  # individual power
  power.ind <- apply(rejects.nonzero, 2, mean)
  power.ind.mean <- mean(power.ind)

  # minimum and complete power
  power.min <- rep(NA, num.nonzero)

  # if unadjusted, don't report minimum or complete power
  if(adj)
  {
    for(m in 1:num.nonzero)
    {
      min.rejects <- apply(rejects.nonzero, 1, function(x){ sum(x) >= m })
      power.min[m] <- mean(min.rejects)
    }
  }


  if(num.nonzero == 0)
  {
    all.power.results <- data.frame('D1indiv' = NA)
  } else
  {
    # combine all power for all definitions
    all.power.results <- data.frame(matrix(c(power.ind, power.ind.mean, power.min), nrow = 1))

    if(num.nonzero > 1)
    {
      colnames(all.power.results) = c(paste0("D", 1:num.nonzero, "indiv"),
                                      "indiv.mean", paste0("min",1:(num.nonzero-1)), "complete")
    } else
    {
      colnames(all.power.results) = c(paste0("D", 1:num.nonzero, "indiv"),
                                      "indiv.mean", "min1")
    }
  }

  return(all.power.results)
}


#' Calculate power using PUMP method
#'
#' This functions calculates power for all definitions of power (individual,
#' d-minimal, complete) for all the different MTPs (none, Bonferroni, Holms,
#' Bejamini-Hocheberg, Westfall-Young Single Step, Westfall-Young Step Down).
#'
#' @param design a single RCT design (see list/naming convention)
#' @param MTP Multiple adjustment procedure of interest. Supported options:
#'   none, Bonferroni, BH, Holm, WY-SS, WY-SD (passed as strings).  Provide list to
#'   automatically re-run for each procedure on the list.
#' @param MDES scalar or vector:  t he MDES values for each outcome.
#' Please provide a scalar, a vector of length M, or vector of values for non-zero outcomes.
#' @param numZero Additional number of outcomes assumed to be zero. Please provide NumZero + length(MDES) = M
#' @param M scalar; the number of hypothesis tests (outcomes), including zero outcomes
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
  design, MTP = NULL, MDES, numZero = NULL,
  M,
  nbar, J = 1, K = 1, Tbar,
  alpha = 0.05,
  numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
  R2.1 = 0, R2.2 = 0, R2.3 = 0,
  ICC.2 = 0, ICC.3 = 0,
  omega.2 = 0, omega.3 = 0,
  rho = NULL, rho.matrix = NULL,
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
      MTP = MTP,
      MDES = MDES, numZero = numZero, M = M, J = J, K = K,
      nbar = nbar, Tbar = Tbar, alpha = alpha,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
      ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
      rho = rho, rho.matrix = rho.matrix, B = B
    )

    params.list <- validate_inputs(design, params.list)

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
  if(is.null(rho.matrix))
  {
    Sigma <- matrix(rho, M, M)
    diag(Sigma) <- 1
  } else
  {
    Sigma <- rho.matrix
  }

  # generate t values and p values under alternative hypothesis using multivariate t-distribution
  rawt.mat <- mvtnorm::rmvt(tnum, sigma = Sigma, df = t.df) + t.shift.mat
  rawp.mat <- pt(-abs(rawt.mat), df = t.df) * 2

  if (is.function(updateProgress) & !is.null(rawp.mat)) {
    updateProgress(message = "P-values have been generated!")
  }

  grab.pval <- function(...,proc) {return(...$adjp[order(...$index),proc])}

  adjp = NULL
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

    adjp <- adjp.wyss(rawt.mat = rawt.mat, B = B, Sigma = Sigma, t.df = t.df)

  } else if (MTP == "WY-SD"){

    adjp <- adjp.wysd(rawt.mat = rawt.mat, B = B, Sigma = Sigma, t.df = t.df, cl = cl)

  } else if ( MTP != "None") {

    stop(paste("Unknown MTP:", MTP))
  }

  if (is.function(updateProgress) & !is.null(adjp)){
    updateProgress(message = paste("Multiple adjustments done for", MTP))
  }

  ind.nonzero <- MDES > 0
  power.results.rawp <- get.power.results(rawp.mat, ind.nonzero, alpha, adj = FALSE)

  if ( !is.null( adjp ) ) {
    power.results.proc <- get.power.results(adjp, ind.nonzero, alpha, adj = TRUE)
    power.results.all <- data.frame(rbind(power.results.rawp, power.results.proc))
    rownames(power.results.all) <- c('rawp', MTP)

    return(power.results.all)
  } else {
    rownames(power.results.rawp) = "rawp"
    return(power.results.rawp)
  }
}

