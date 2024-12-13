
#' Calculates p-values from t-values
#'
#' @param rawt vector of t statistics
#' @param t.df degrees of freedom
#' @param two.tailed whether to calculate 1 or 2-tailed p-values
#'
#' @return power results for individual, minimum, complete power
#' @keywords internal
calc_pval <- function(rawt, t.df, two.tailed)
{
    if (two.tailed)
    {
        rawp <- 2*(1 - stats::pt(abs(rawt), df = t.df))
    } else
    {
        rawp <- 1 - stats::pt(rawt, df = t.df)
    }
    return(rawp)
}



#' @title Calculates different definitions of power (support function)
#'
#' @description This function takes in a matrix of adjusted p-values
#'   and unadjusted p-values and outputs different types of power.
#'
#'   This function is mostly for internal use, but may be of interest
#'   to users who wish to calculate power on their own.
#'
#' @param adj.pval.mat matrix; adjusted p-values, columns are outcomes
#' @param unadj.pval.mat matrix; unadjusted p-values, columns are
#'   outcomes
#' @param ind.nonzero vector; which outcomes correspond to nonzero
#'   effects.
#' @param alpha scalar; the family wise error rate (FWER).
#' @param drop.zero.outcomes logical; whether to report power results
#'   for outcomes with MDES = 0.
#' @param adj logical; whether p-values are unadjusted or not.
#'
#' @return data frame; power results for individual, minimum, complete
#'   power.
#'
#' @export
get_power_results <- function(adj.pval.mat, unadj.pval.mat,
                              ind.nonzero, alpha,
                              drop.zero.outcomes = TRUE,
                              adj = TRUE)
{
    M <- ncol(adj.pval.mat)
    num.nonzero <- sum(ind.nonzero)
    
    # if all zeros, then return first thing.
    if (num.nonzero == 0 && drop.zero.outcomes) {
        #all.power.results <- data.frame('D1indiv' = NA)
        ind.nonzero[1] = TRUE
    }
    
    # rejected tests
    rejects <- apply(adj.pval.mat, 2, function(x){ 1*(x < alpha) })
    # unadjusted
    rejects.unadj <- apply(unadj.pval.mat, 2, function(x){ 1*(x < alpha) })
    
    # individual power
    power.ind <- apply(rejects, 2, mean)
    names(power.ind) <- paste0('D', 1:M, 'indiv')
    
    # minimum power
    power.min <- rep(NA, M)
    names(power.min) <- paste0('min',1:M)
    
    # complete power
    power.complete <- NA
    
    # if unadjusted, don't report minimum or complete power
    if (adj)
    {
        # complete power
        if (all(ind.nonzero))
        {
            complete.rejects <- apply(rejects.unadj, 1, 
                                      function(x){ sum(x) == M })
            power.complete <- mean(complete.rejects)
        }
        # minimum power
        for (m in 1:M)
        {
            min.rejects <- apply(rejects, 1, function(x){ sum(x) >= m })
            power.min[m] <- mean(min.rejects)
        }
    }
    
    # subset to only nonzero where relevant
    if (drop.zero.outcomes)
    {
        power.ind <- power.ind[ind.nonzero]
        power.min <- power.min[ind.nonzero]
        # rename to be more sensible if there are zeros in the middle
        names(power.min) <- paste0('min',seq(1, length(power.min)))
    }
    
    power.ind.mean <- mean(power.ind)
    names(power.ind.mean) <- 'indiv.mean'
    names(power.complete) <- 'complete'
    
    # remove redundant min column
    if (sum(num.nonzero) == M)
    {
        power.min <- power.min[1:(M - 1)]
    }
    
    if (M > 1)
    {
        power.vec <- c(power.ind, power.ind.mean, power.min, power.complete)
        # don't return min and complete for M = 1
    } else
    {
        power.vec <- c(power.ind)
    }
    
    power.vec <- vapply(power.vec, as.numeric, numeric(1))
    
    # combine all power for all definitions
    all.power.results <- data.frame(matrix(power.vec, nrow = 1))
    colnames(all.power.results) <- names(power.vec)
    
    
    return(all.power.results)
}




calc_MT <- function(df, alpha, two.tailed, target.power) {
    # t statistics
    
    T2 <- abs(stats::qt(target.power, df))
    
    # number of SEs we need for MDES
    MT <- ifelse(target.power >= 0.5, T1 + T2, T1 - T2)
    
    return(MT)
}


#' Calculate power theoretically for M=1 situations
#'
#' @param MDES MDES (single number)
#' @param SE Calculated SE of the estimator
#' @param df Degrees of freedom of the estimator (often approximated).
#' @param alpha Alpha for the planned test.
#' @param two.tailed TRUE/FALSE  Two- or one-sided test?
#' @return Single row Tibble with columns or power, SE, and DF.  MTP
#'   column with value of "None".
#'
#' @keywords internal
pump_power_exact <- function( MDES, SE, df, alpha, two.tailed ) {
    
    T1 <- ifelse(two.tailed == TRUE, 
                 abs(stats::qt(alpha/2, df)), 
                 abs(stats::qt(alpha, df)))
    
    pow <- stats::pt( T1 - MDES / SE, df = df, lower.tail = FALSE ) +
        stats::pt( -T1 - MDES / SE, df = df, lower.tail = TRUE )
    
    dplyr::tibble( MTP = "None", D1indiv = pow, SE1 = SE, df1 = df )
}



#' @title Estimate power across definitions (core function)
#'
#' @description The user chooses the context (d_m), MTP, MDES, and
#'   choices of all relevant design parameters.
#'
#'   The functions returns power for all definitions of power for any
#'   MTP. For a list of choices for specific parameters, see
#'   pump_info().
#'
#' @seealso For more detailed information about this function and the
#'   user choices, see the manuscript <doi:10.18637/jss.v108.i06>,
#'   which includes a detailed Technical Appendix including
#'   information about the designs and models and parameters.
#'
#' @param d_m string; a single context, which is a design and model
#'   code. See pump_info() for list of choices.
#' @param MTP string, or vector of strings; multiple testing
#'   procedure(s). See pump_info() for list of choices.
#' @param MDES scalar or vector; the desired MDES values for each
#'   outcome. Please provide a scalar, a vector of length M, or vector
#'   of values for non-zero outcomes.
#' @param numZero scalar; additional number of outcomes assumed to be
#'   zero. Please provide NumZero + length(MDES) = M, if length(MDES)
#'   is not 1.
#' @param propZero scalar; proportion of outcomes assumed to be zero
#'   (alternative specification to numZero).  length(MDES) should be 1
#'   or equal to (1-propZero)*M.
#' @param M scalar; the number of hypothesis tests (outcomes),
#'   including zero outcomes.
#' @param J scalar; the harmonic mean of number of level 2 units per
#'   level 3 unit (schools per district). Note that this is not the
#'   total number of level 2 units, but instead the number of level 2
#'   units nested within each level 3 unit, so the total number of
#'   level 2 units is J x K.
#' @param K scalar; the number of level 3 units (districts).
#' @param nbar scalar; the harmonic mean of the number of level 1
#'   units per level 2 unit (students per school). Note that this is
#'   not the total number of level 1 units, but instead the number of
#'   level 1 units nested within each level 2 unit, so the total
#'   number of level 1 units is nbar x J x K.
#' @param Tbar scalar; the proportion of samples that are assigned to
#'   the treatment.
#' @param alpha scalar; the family wise error rate (FWER).
#' @param two.tailed scalar; TRUE/FALSE for two-tailed or one-tailed
#'   power calculation.
#' @param numCovar.1 scalar; number of level 1 (individual)
#'   covariates.
#' @param numCovar.2 scalar; number of level 2 (school) covariates.
#' @param numCovar.3 scalar; number of level 3 (district) covariates.
#' @param R2.1 scalar, or vector of length M; percent of variation
#'   explained by level 1 covariates for each outcome.
#' @param R2.2 scalar, or vector of length M; percent of variation
#'   explained by level 2 covariates for each outcome.
#' @param R2.3 scalar, or vector of length M; percent of variation
#'   explained by level 3 covariates for each outcome.
#' @param ICC.2 scalar, or vector of length M; level 2 (school)
#'   intraclass correlation.
#' @param ICC.3 scalar, or vector length M; level 3 (district)
#'   intraclass correlation.
#' @param omega.2 scalar, or vector of length M; ratio of variance of
#'   level 2 average impacts to variance of level 2 random intercepts.
#' @param omega.3 scalar, or vector of length M; ratio of variance of
#'   level 3 average impacts to variance of level 3 random intercepts.
#' @param rho scalar; assumed correlation between all pairs of test
#'   statistics.
#' @param rho.matrix matrix; alternate specification allowing a full
#'   matrix of correlations between test statistics. Must specify
#'   either rho or rho.matrix, but not both.
#' @param tnum scalar; the number of test statistics to draw.
#'   Increasing tnum increases precision and computation time.
#' @param B scalar; the number of permutations for Westfall-Young
#'   procedures.
#' @param parallel.WY.cores number of cores to use for parallel
#'   processing of WY-SD.
#' @param drop.zero.outcomes whether to report power results for
#'   outcomes with MDES = 0.  If ALL MDES = 0, then the first outcome
#'   will not be dropped.
#' @param updateProgress function to update progress bar (only used
#'   for PUMP shiny app).
#' @param long.table TRUE for table with power as rows, correction as
#'   columns, and with more verbose names. See
#'   `transpose_power_table`.
#' @param verbose TRUE/FALSE; Print out diagnostics of time, etc.
#' @param validate.inputs TRUE/FALSE; whether or not to check whether
#'   parameters are valid given the choice of d_m.
#' @param exact.where.possible TRUE/FALSE; whether to do exact
#'   calculations when M=1, or use simulation.  Default is TRUE.
#' @return a pumpresult object containing power results.
#' @export
#'
#' @examples
#' pp <- pump_power(
#'    d_m = "d3.2_m3ff2rc",
#'    MTP = 'HO',
#'    nbar = 50,
#'    J = 30,
#'    K = 10,
#'    M = 5,
#'    MDES = 0.125,
#'    Tbar = 0.5, alpha = 0.05,
#'    numCovar.1 = 1, numCovar.2 = 1,
#'    R2.1 = 0.1, R2.2 = 0.1,
#'    ICC.2 = 0.2, ICC.3 = 0.2,
#'    omega.2 = 0, omega.3 = 0.1,
#'    rho = 0.5)
#' 
pump_power <- function(
        d_m, MTP = NULL, MDES, 
        numZero = NULL, propZero = NULL,
        M = 1,
        nbar, J = 1, K = 1, Tbar,
        alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
        R2.1 = 0, R2.2 = 0, R2.3 = 0,
        ICC.2 = 0, ICC.3 = 0,
        omega.2 = 0, omega.3 = 0,
        rho = NULL, rho.matrix = NULL,
        tnum = 10000, B = 1000,
        parallel.WY.cores = 1,
        drop.zero.outcomes = TRUE,
        updateProgress = NULL,
        validate.inputs = TRUE,
        long.table = FALSE,
        verbose = FALSE,
        exact.where.possible = TRUE
) {
    # do not duplicate 'None'
    if (length(MTP) > 1 && any(MTP == 'None')) {
        MTP <- MTP[which(MTP != 'None')]
    }
    
    
    # Recursively call self for each element on MTP list.
    MTP_list = NULL
    if ( length( MTP ) > 1 ) {
        if ( verbose ) {
            smessage( "Multiple MTPs leading to %d calls\n", length(MTP) )
        }
        MTP_list <- validate_MTP( MTP = MTP,
                      power.call = TRUE, mdes.call = FALSE, ss.call = FALSE,
                      M = M, pdef = NULL, multi.MTP.ok = TRUE )
        
    }
    
    
    if (validate.inputs) {
        # validate input parameters
        params.list <- list(
            MTP = MTP[1],
            MDES = MDES, numZero = numZero, propZero = propZero, 
            M = M, J = J, K = K,
            nbar = nbar, Tbar = Tbar, alpha = alpha, two.tailed = two.tailed,
            numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, 
            numCovar.3 = numCovar.3,
            R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
            ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
            rho = rho, rho.matrix = rho.matrix, B = B, tnum = tnum,
            power.definition = NULL
        )
        params.list <- validate_inputs(
            d_m, params.list, power.call = TRUE, verbose = verbose 
        )
        
        MTP <- params.list$MTP
        MDES <- params.list$MDES
        M <- params.list$M; J <- params.list$J; K <- params.list$K
        nbar <- params.list$nbar; Tbar <- params.list$Tbar
        alpha <- params.list$alpha; two.tailed <- params.list$two.tailed
        numCovar.1 <- params.list$numCovar.1; numCovar.2 <- params.list$numCovar.2
        numCovar.3 <- params.list$numCovar.3
        R2.1 <- params.list$R2.1; R2.2 <- params.list$R2.2
        R2.3 <- params.list$R2.3
        ICC.2 <- params.list$ICC.2; ICC.3 <- params.list$ICC.3
        omega.2 <- params.list$omega.2; omega.3 <- params.list$omega.3
        rho <- params.list$rho; rho.matrix <- params.list$rho.matrix
        d_m = params.list$d_m
        
        params.list <- params.list[names(params.list) != 'power.definition']
    } else {
        params.list <- NULL
    }
    
    if ( !is.null( MTP_list ) ) {
        des <- purrr::map( MTP_list,
                            pump_power, 
                            d_m = d_m, MDES = MDES,
                            M = M, J = J, K = K, nbar = nbar,
                            numZero = numZero, propZero = propZero,
                            Tbar = Tbar,
                            alpha = alpha, two.tailed = two.tailed,
                            numCovar.1 = numCovar.1, numCovar.2 = numCovar.2,
                            numCovar.3 = numCovar.3,
                            R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                            ICC.2 = ICC.2, ICC.3 = ICC.3,
                            rho = rho, rho.matrix = rho.matrix,
                            omega.2 = omega.2, omega.3 = omega.3,
                            long.table = long.table,
                            tnum = tnum, B = B, 
                            parallel.WY.cores = parallel.WY.cores,
                            drop.zero.outcomes = drop.zero.outcomes,
                            verbose = verbose,
                            exact.where.possible = exact.where.possible,
                            validate.inputs = FALSE,
                            updateProgress = updateProgress )
        
        plist <- params.list
        plist$MTP <- MTP_list
        
        if ( long.table ) {
            ftable <- des[[1]]
            for (i in 2:length(des)) {
                ftable <- dplyr::bind_cols( ftable, des[[i]][ ncol(des[[i]]) ] )
            }
        } else {
            ftable <- des[[1]]
            for (i in 2:length(des)) {
                ftable <- dplyr::bind_rows( ftable, des[[i]][ nrow(des[[i]]), ] )
            }
        }
        
        return( make.pumpresult( ftable, "power",
                                 params.list = plist,
                                 d_m = d_m,
                                 long.table = long.table ) )
    }
    
    
    
    
    # compute test statistics for when null hypothesis is false
    Q.m <- calc_SE(
        d_m = d_m, J = J, K = K, nbar = nbar, Tbar = Tbar,
        R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
        ICC.2 = ICC.2, ICC.3 = ICC.3,
        omega.2 = omega.2, omega.3 = omega.3
    )
    t.shift <- MDES/Q.m
    t.df <- calc_df(
        d_m = d_m, J = J, K = K,
        nbar = nbar,
        numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
        validate = validate.inputs
    )
    t.shift.mat <- t(matrix(rep(t.shift, tnum), M, tnum))
    
    if ( M == 1 && exact.where.possible ) {
        ppe <- pump_power_exact( MDES, Q.m, t.df,
                                 alpha = alpha,
                                 two.tailed = two.tailed )
        
        if ( long.table ) {
            power.results <- transpose_power_table( ppe, M = M )
        }
        return( make.pumpresult( ppe,
                                 "power",
                                 params.list = params.list,
                                 d_m = d_m,
                                 exact = TRUE,
                                 long.table = long.table ) )
    }
    
    
    # correlation between the test statistics
    if (is.null(rho.matrix))
    {
        Sigma <- matrix(rho, M, M)
        diag(Sigma) <- 1
    } else
    {
        Sigma <- rho.matrix
    }
    
    
    # generate t values and p values under alternative hypothesis
    # using multivariate t-distribution
    rawt.mat <- matrix(mvtnorm::rmvt(tnum, sigma = Sigma, df = t.df) + 
                           t.shift.mat, nrow = tnum, ncol = M)
    rawp.mat <- calc_pval(rawt.mat, t.df, two.tailed)
    
    if (is.function(updateProgress) & !is.null(rawp.mat)) {
        updateProgress(message = "P-values have been generated!")
    }
    
    if (MTP == "BF") {
        
        adjp.mat <- t(apply(rawp.mat, 1, stats::p.adjust, method = "bonferroni"))
        
    } else if (MTP == "HO") {
        
        adjp.mat <- t(apply(rawp.mat, 1, stats::p.adjust, method = "holm"))
        
    } else if (MTP == "BH") {
        
        adjp.mat <- t(apply(rawp.mat, 1, stats::p.adjust, method = "BH"))
        
    } else if (MTP == "WY-SS") {
        
        adjp.mat <- adjp_wyss(rawp.mat = rawp.mat, B = B,
                              Sigma = Sigma, t.df = t.df,
                              two.tailed = two.tailed,
                              updateProgress = updateProgress)
        
    } else if (MTP == "WY-SD") {
        
        if ( parallel.WY.cores > 1 )
        {
            cl <- parallel::makeCluster(parallel.WY.cores)
        } else
        {
            cl <- NULL
        }
        
        adjp.mat <- adjp_wysd(rawp.mat = rawp.mat, B = B,
                              Sigma = Sigma, t.df = t.df,
                              two.tailed = two.tailed, cl = cl,
                              updateProgress = updateProgress)
        
        if ( parallel.WY.cores > 1 )
        {
            parallel::stopCluster(cl)
        }
        
    } else {
        adjp.mat <- NULL
    }
    
    if (is.function(updateProgress) & !is.null(adjp.mat)) {
        updateProgress(message = paste("Multiple adjustments done for", MTP))
    }
    
    ind.nonzero <- MDES > 0
    
    power.results.raw <- get_power_results(
        adj.pval.mat = rawp.mat, unadj.pval.mat = rawp.mat,
        ind.nonzero = ind.nonzero, alpha = alpha,
        drop.zero.outcomes = drop.zero.outcomes, adj = FALSE
    )
    
    
    SEs <- data.frame(t(c(Q.m, t.df)))
    colnames(SEs) <- c(paste0("SE", 1:length(Q.m)), paste0("df", 1:length(t.df)))
    
    power.results.raw = dplyr::bind_cols(power.results.raw, SEs)
    
    if ( MTP != 'None' ) {
        power.results.proc <- get_power_results(
            adj.pval.mat = adjp.mat, unadj.pval.mat = rawp.mat,
            ind.nonzero = ind.nonzero, alpha = alpha,
            drop.zero.outcomes = drop.zero.outcomes, adj = TRUE
        )
        
        power.results <- dplyr::bind_rows( None = power.results.raw,
                                           Alt = power.results.proc,
                                           .id = "MTP" )
        power.results$MTP[[2]] = MTP
    } else {
        power.results <- cbind('MTP' = 'None', power.results.raw)
    }
    
    
    if ( long.table ) {
        power.results <- transpose_power_table( power.results, M = M )
    }
    return( make.pumpresult( power.results,
                             "power",
                             params.list = params.list,
                             d_m = d_m,
                             long.table = long.table ) )
}
