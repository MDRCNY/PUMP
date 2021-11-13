#' Convert power table from wide to long
#'
#' Transform table returned from pump_power to a long format table.
#'
#' @param power_table pumpresult object for a power result (not mdes or sample).
#' 
transpose_power_table <- function( power_table ) {

  cname = power_table$MTP
  power_table$MTP = NULL
  pp <- t( power_table )
  colnames(pp) <- cname
  #if ( ncol( pp ) > 1 ) {
  #  pp = pp[ , ncol(pp):1 ]
  #}
  pows <- rownames(pp)
  pp <- pp %>% # pp[ nrow(pp):1, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column( var="power" )

  pp$power <- stringr::str_replace( pp$power, "D(.*)indiv", "individual outcome \\1" )
  pp$power <- stringr::str_replace( pp$power, "min(.*)", "\\1-minimum" )
  pp$power <- stringr::str_replace( pp$power, "indiv.mean", "mean individual" )
  pp
}


#' Calculates different definitions of power
#'
#' This function takes in a matrix of adjusted p-values and outputs different types of power
#'
#' @param pval.mat matrix of p-values, columns are outcomes
#' @param ind.nonzero vector indicating which outcomes are nonzero
#' @param alpha scalar; the family wise error rate (FWER)
#' @param adj whether p-values are unadjusted or not
#'
#' @return power results for individual, minimum, complete power
#' @export
get.power.results <- function(pval.mat, ind.nonzero, alpha, adj = TRUE)
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
#'   none, Bonferroni, BH, Holm, WY-SS, WY-SD (passed as strings).  Provide list
#'   to automatically re-run for each procedure on the list.
#' @param MDES scalar or vector; the desired MDES values for each outcome. Please
#'   provide a scalar, a vector of length M, or vector of values for non-zero
#'   outcomes.
#' @param numZero scalar; Additional number of outcomes assumed to be zero. Please
#'   provide NumZero + length(MDES) = M
#' @param M scalar; the number of hypothesis tests (outcomes), including zero
#'   outcomes
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
#' @param ICC.2 scalar, or vector of length M; school intraclass correlation
#' @param ICC.3 scalar, or vector length M; district intraclass correlation
#' @param omega.2 scalar, or vector of length M; ratio of variance of school-average impacts to
#'   variance of school-level random intercepts.  Default to 0 (no treatment
#'   variation).
#' @param omega.3 scalar, or vector of length M; ratio of variance of district-average impacts to
#'   variance of district-level random intercepts. Default to 0 (no treatment
#'   variation).
#' @param rho scalar; assumed correlation between all pairs of test statistics.
#' @param rho.matrix matrix; alternate specification allowing a full matrix
#' of correlations between test statistics. Must specify either
#' rho or rho.matrix, but not both.
#' @param tnum scalar; the number of test statistics (samples)
#' @param B scalar; the number of samples/permutations for Westfall-Young
#' @param cl cluster object to use for parallel processing
#' @param updateProgress the callback function to update the progress bar (User
#'   does not have to input anything)
#' @param long.table TRUE for table with power as rows, correction as columns,
#'   and with more verbose names.  See `transpose_power_table`.
#' @param verbose Print out diagnostics of time, etc.
#' @param validate.inputs whether or not to check whether parameters are valid
#' given the choice of design
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
  validate.inputs = TRUE,
  long.table = FALSE,
  verbose = FALSE
)
{
  # Call self for each element on MTP list.
  if ( length( MTP ) > 1 ) {
    if ( verbose ) {
      scat( "Multiple MTPs leading to %d calls\n", length(MTP) )
    }
    des = purrr::map( MTP,
                      pump_power, design = design, MDES = MDES,
                      M = M, J = J, K = K, nbar = nbar, numZero = numZero,
                      Tbar = Tbar,
                      alpha = alpha,
                      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2,
                      numCovar.3 = numCovar.3,
                      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                      ICC.2 = ICC.2, ICC.3 = ICC.3,
                      rho = rho, omega.2 = omega.2, omega.3 = omega.3,
                      long.table = long.table,
                      tnum = tnum, B = B, cl = cl,
                      updateProgress = updateProgress )

    plist = attr( des[[1]], "params.list" )
    plist$MTP = MTP
    if ( long.table ) {
      ftable = des[[1]]
      for ( i in 2:length(des) ) {
        ftable = dplyr::bind_cols( ftable, des[[i]][ ncol(des[[i]]) ] )
      }
    } else {
      ftable = des[[1]]
      for ( i in 2:length(des) ) {
        ftable = dplyr::bind_rows( ftable, des[[i]][ nrow(des[[i]]), ] )
      }
    }
    return( make.pumpresult( ftable, "power",
                             params.list = plist,
                             design = design,
                             multiple_MTP = TRUE,
                             long.table=long.table ) )

    #des = map( des, ~ .x[nrow(.x),] ) %>%
    #  dplyr::bind_rows()
    #return( des )
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

    params.list <- validate_inputs(design, params.list, power.call = TRUE)

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
  } else {
    params.list = NULL
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
  rawp.mat <- stats::pt(-abs(rawt.mat), df = t.df) * 2

  if (is.function(updateProgress) & !is.null(rawp.mat)) {
    updateProgress(message = "P-values have been generated!")
  }

  if (MTP == "Bonferroni"){

    adjp <- t(apply(rawp.mat, 1, stats::p.adjust, method = "bonferroni"))

  } else if (MTP == "Holm") {

    adjp <- t(apply(rawp.mat, 1, stats::p.adjust, method = "holm"))

  } else if (MTP == "BH") {

    adjp <- t(apply(rawp.mat, 1, stats::p.adjust, method = "hochberg"))

  } else if (MTP == "WY-SS"){

    adjp <- adjp.wyss(rawt.mat = rawt.mat, B = B,
                      Sigma = Sigma, t.df = t.df)

  } else if (MTP == "WY-SD"){

    adjp <- adjp.wysd(rawt.mat = rawt.mat, B = B,
                      Sigma = Sigma, t.df = t.df, cl = cl)

  } else
  {
    adjp <- NULL
  }

  if (is.function(updateProgress) & !is.null(adjp)){
    updateProgress(message = paste("Multiple adjustments done for", MTP))
  }

  ind.nonzero <- MDES > 0
  power.results.raw <- get.power.results(rawp.mat, ind.nonzero, alpha, adj = FALSE)

  if ( MTP != 'None' ) {
    power.results.proc <- get.power.results(adjp, ind.nonzero, alpha, adj = TRUE)
    power.results <- data.frame(rbind(power.results.raw, power.results.proc))
    power.results <- cbind('MTP' = c('None', MTP), power.results)
  } else {
    power.results <- cbind('MTP' = 'None', power.results.raw)
  }
  if ( long.table ) {
    power.results <- transpose_power_table( power.results )
  }
  return( make.pumpresult( power.results, "power",
                           params.list = params.list,
                           design = design,
                           long.table = long.table ) )
}
