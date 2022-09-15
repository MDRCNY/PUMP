
# Code for the pump_sample method

calc_MT <- function( df, alpha, two.tailed, target.power ) {
    # t statistics
    T1 <- ifelse(two.tailed == TRUE, 
                 abs(stats::qt(alpha/2, df)), 
                 abs(stats::qt(alpha, df)))
    T2 <- abs(stats::qt(target.power, df))
    
    # number of SEs we need for MDES
    MT <- ifelse(target.power >= 0.5, T1 + T2, T1 - T2)
    
    return(MT)
}



#' Calculating Needed Sample Size for Raw (Unadjusted) Power
#'
#' This is a helper function for getting a needed Sample Size when no
#' adjustments has been made to the test statistics.
#'
#' It is for a single, individual outcome.  It only takes scalar values for all
#' its arguments, and does not have an M argument (for number of outcomes).
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
#' @param typesample type of sample size to calculate: J, K, or nbar
#' @param target.power target power to arrive at
#' @param max.steps how many steps allowed before terminating
#' @param warn.small Warn if degrees of freedom issues are causing inability to
#'   achieve target power for sample size.
#'
#' @return Requisit sample size (as integer) and associated degreess of freedom.
#' @keywords internal
pump_sample_raw <- function(
        d_m, MTP, typesample,
        MDES,
        nbar = NULL, J = NULL, K = NULL,
        target.power,
        Tbar, alpha = 0.05, two.tailed,
        numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
        R2.1, R2.2 = NULL, R2.3 = NULL, 
        ICC.2 = NULL, ICC.3 = NULL,
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
    
    
    # check for vectorized components
    if(length(MDES) > 1 |
       length(R2.1) > 1 | length(R2.2) > 1 | length(R2.3) > 1 |
       length(ICC.2) > 1 | length(ICC.3) > 1 |
       length(omega.2) > 1 | length(omega.3) > 1
    )
    {
        stop('pump_sample_raw only takes scalar inputs')
    }
    
    initial_df <- calc_df(d_m, J, K, nbar, 
                          numCovar.1, numCovar.2, numCovar.3)
    
    stopifnot( initial_df > 0 )
    
    i <- 0
    conv <- FALSE
    
    # Get initial size (will be low)
    MT <- calc_MT(df = initial_df, alpha = alpha, 
                  two.tailed = two.tailed, 
                  target.power = target.power)
    if (typesample == "J") {
        J <- calc_J( d_m, MT = MT, MDES = MDES[1],
                     K = K, nbar = nbar, Tbar = Tbar,
                     R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
                     ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
                     omega.2 = omega.2[1], omega.3 = omega.3[1] )
        J <- round(J)
    } else if (typesample == "K") {
        K <- calc_K(
            d_m, MT = MT, MDES = MDES[1],
            J = J, nbar = nbar, Tbar = Tbar,
            R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
            ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
            omega.2 = omega.2[1], omega.3 = omega.3[1]
        )
        K <- round(K)
    } else if (typesample == "nbar") {
        nbar <- calc_nbar(
            d_m, MT = MT, MDES = MDES[1], 
            J = J, K = K, Tbar = Tbar,
            R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
            ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
            omega.2 = omega.2[1], omega.3 = omega.3[1]
        )
    }
    
    df <- calc_df(d_m, J, K, nbar, 
                  numCovar.1, numCovar.2, numCovar.3, 
                  validate = FALSE)
    
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
            
            df <- calc_df(d_m, J, K, nbar, 
                          numCovar.1, numCovar.2, numCovar.3, 
                          validate = FALSE)
            
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
        
        df <- calc_df(d_m, J, K, nbar, 
                      numCovar.1, numCovar.2, numCovar.3, 
                      validate = FALSE)
        MT <- calc_MT(df = df, alpha = alpha, 
                      two.tailed = two.tailed, 
                      target.power = target.power)
        
        
        if (typesample == "J") {
            J1 <- calc_J( d_m, MT = MT, MDES = MDES[1],
                          K = K, nbar = nbar, Tbar = Tbar,
                          R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
                          ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
                          omega.2 = omega.2[1], omega.3 = omega.3[1] )
            J1 <- round( J1 )
            
            if ( is.na(J1) || (J1 <= J) ) {
                conv <- TRUE
            } else {
                J <- J + 1
            }
            
        } else if (typesample == "K") {
            K1 <- calc_K(
                d_m, MT = MT, MDES = MDES[1],
                J = J, nbar = nbar, Tbar = Tbar,
                R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
                ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
                omega.2 = omega.2[1], omega.3 = omega.3[1]
            )
            K1 <- round( K1 )
            if ( is.na(K1) || (K1 <= K) ) {
                conv <- TRUE
            } else {
                K <- K + 1
            }
        } else if (typesample == "nbar") {
            nbar1 <- calc_nbar(
                d_m, MT = MT, MDES = MDES[1],
                J = J, K = K, Tbar = Tbar,
                R2.1 = R2.1[1], R2.2 = R2.2[1], R2.3 = R2.3[1],
                ICC.2 = ICC.2[1], ICC.3 = ICC.3[1],
                omega.2 = omega.2[1], omega.3 = omega.3[1]
            )
            
            if (is.na( nbar1 ) || (nbar1 <= nbar) ) {
                conv <- TRUE
            } else {
                nbar <- nbar + 1
            }
        }
        
        i <- i + 1
    }
    
    if ( i >= max.steps ) {
        stop( "Hit maximum iterations in pump_sample_raw()" )
    }
    
    if (typesample == "J") {
        if(!is.na(J) & J <= 0){ J <- NA }
        return( list( ss = J, df = df ) )
    } else if (typesample == "K") {
        if(!is.na(K) & K <= 0){ K <- NA }
        return( list( ss = K, df = df ) )
    } else if (typesample == "nbar") {
        if(!is.na(nbar) & nbar <= 0){ nbar <- NA }
        return( list( ss = nbar, df = df ) )
    }
}






#' @title Estimate the required sample size (core function)
#'
#' @description The user chooses the context (d_m), MTP,
#' type of sample size, 
#' MDES,
#' power definition, and choices of all relevant design parameters.
#' 
#' The functions performs a search algorithm,
#' and returns the sample size value within the specified tolerance.
#' For a list of choices for specific parameters, see pump_info().
#' 
#' @seealso For more detailed information about this function 
#' and the user choices,
#' see the manuscript \url{https://arxiv.org/abs/2112.15273},
#' which includes a detailed Technical Appendix
#' including information about the designs and models
#' and parameters.
#'
#' @inheritParams pump_mdes
#' @inheritParams pump_power
#'
#' @param typesample string; type of sample size to 
#' calculate: "nbar", "J", or "K".
#' @param max_sample_size_nbar scalar; default upper bound for nbar 
#' for search algorithm.
#' @param max_sample_size_JK scalar; default upper bound for J or K 
#' for search algorithm.
#'
#' @return a pumpresult object containing sample size results.
#' @export
#' 
#' @examples
#' J <- pump_sample(
#'   d_m = 'd2.1_m2fc',
#'   MTP = 'HO',
#'   power.definition = 'D1indiv',
#'   typesample = 'J',
#'   target.power = 0.8,
#'   nbar = 50,
#'   M = 3,
#'   MDES = 0.125,
#'   Tbar = 0.5, alpha = 0.05,
#'   numCovar.1 = 1,
#'   R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2,
#'   tnum = 1000)

pump_sample <- function(
        d_m, MTP = NULL, typesample,
        MDES, M = 1, numZero = NULL,
        nbar = NULL, J = NULL, K = NULL,
        target.power, power.definition,
        alpha, two.tailed = TRUE,
        Tbar,
        numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
        R2.1 = 0, R2.2 = 0, R2.3 = 0,
        ICC.2 = 0, ICC.3 = 0,
        rho = NULL, rho.matrix = NULL,
        omega.2 = 0, omega.3 = 0,
        B = 1000,
        max.steps = 20, tnum = 1000, start.tnum = tnum / 10, final.tnum = 4*tnum,
        parallel.WY.cores = 1, updateProgress = NULL,
        max_sample_size_nbar = 10000,
        max_sample_size_JK = 1000,
        tol = 0.01, give.optimizer.warnings = FALSE,
        verbose = FALSE )
{
    if ( M == 1 ) {
        if ( missing( "power.definition" ) ) {
            power.definition = "D1indiv"
        } 
    }
    
    if ( verbose ) {
        scat( "pump_mdes with %d max iterations per search, 
          starting at %d iterations with final %d iterations.
          \n\tMax steps %d\n\t%d perms for WY if used\n",
              tnum, start.tnum, final.tnum, max.steps, B )
    }
    
    
    # Give prelim values for the validation of parameters process.
    if ( typesample == "nbar" ) {
        if(!is.null(nbar)) {
            stop('Do not provide nbar if you are searching for nbar')
        }
        nbar <- 1000
    } else if ( typesample == "J" ) {
        if(!is.null(J)) {
            stop('Do not provide J if you are searching for J')
        }
        J <- 1000
    } else if ( typesample == "K" ) {
        if(!is.null(K)) {
            stop('Do not provide K if you are searching for K')
        }
        K <- 1000
    } else {
        stop( glue::glue( "Invalid typesample '{typesample}'" ) )
    }
    
    # validate input parameters
    params.list <- list(
        MTP = MTP, MDES = MDES, M = M, J = J, K = K, numZero = numZero, 
        nbar = nbar, Tbar = Tbar, alpha = alpha, two.tailed = two.tailed,
        numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
        R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
        ICC.2 = ICC.2, ICC.3 = ICC.3, omega.2 = omega.2, omega.3 = omega.3,
        rho = rho, rho.matrix = rho.matrix, B = B,
        max.steps = max.steps, 
        start.tnum = start.tnum, tnum = tnum, final.tnum = final.tnum,
        power.definition = power.definition
    )
    ##
    params.list <- validate_inputs(
        d_m, params.list, ss.call = TRUE, verbose = verbose 
    )
    ##
    MTP <- params.list$MTP
    MDES <- params.list$MDES; numZero <- params.list$numZero
    M <- params.list$M; J <- params.list$J; K <- params.list$K
    nbar <- params.list$nbar; Tbar <- params.list$Tbar
    alpha <- params.list$alpha; two.tailed <- params.list$two.tailed
    numCovar.1 <- params.list$numCovar.1; numCovar.2 <- params.list$numCovar.2
    numCovar.3 <- params.list$numCovar.3
    R2.1 <- params.list$R2.1; R2.2 <- params.list$R2.2; R2.3 <- params.list$R2.3
    ICC.2 <- params.list$ICC.2; ICC.3 <- params.list$ICC.3
    omega.2 <- params.list$omega.2; omega.3 <- params.list$omega.3
    rho <- params.list$rho; rho.matrix <- params.list$rho.matrix
    B <- params.list$B
    power.definition <- params.list$power.definition
    params.list <- params.list[names(params.list) != 'power.definition']
    
    if ( is.null( numZero ) ) {
        numZero <- 0
        stopifnot( M > numZero )
    }
    
    # power definition type
    pdef <- parse_power_definition( power.definition, M )
    
    pow_params <- list( target.power = target.power,
                        power.definition = power.definition,
                        tol = tol )
    
    # Delete parameter we are actually going to search over.
    if ( typesample == "nbar" ) {
        nbar <- NULL
        params.list["nbar"] <- NULL
    } else if ( typesample == "J" ) {
        J <- NULL
        params.list["J"] <- NULL
    } else if ( typesample == "K" ) {
        K <- NULL
        params.list["K"] <- NULL
    }
    
    output.colnames <- c("MTP", "Sample.type", "Sample.size",
                         paste(power.definition, "power") )
    
    # power checks
    if(round(target.power, 2) <= 0)
    {
        message('Target power of 0 requested')
        ss.results <- data.frame(MTP, typesample, 0, 0)
        colnames(ss.results) <- output.colnames
        return( make.pumpresult( ss.results, type = "sample",
                                 params.list = params.list,
                                 tries = NULL,
                                 d_m = d_m,
                                 sample.level = typesample,
                                 power.params.list = pow_params) )
    }
    if(target.power > 1)
    {
        message('Target power of >1 requested')
        ss.results <- data.frame(MTP, typesample, Inf, 1)
        colnames(ss.results) <- output.colnames
        return( make.pumpresult( ss.results, type = "sample",
                                 params.list = params.list,
                                 tries = NULL,
                                 d_m = d_m,
                                 sample.level = typesample,
                                 power.params.list = pow_params) )
    }
    
    msg <- paste("Estimating sample size of type", typesample, "for",
                 MTP, "for target",
                 power.definition, "power of", round(target.power, 4))
    
    # Checks on what we are estimating, sample size
    if ( verbose ) {
        message(msg)
    }
    if (is.function(updateProgress)) {
        updateProgress(message = msg)
    }
    
    # adjust bounds to capture needed range
    # for minimum or complete power, expand bounds
    # note: complete power is a special case of minimum power
    if ( !pdef$min ) {
        # Compute needed sample size for raw and BF SS for INDIVIDUAL POWER. 
        # We are estimating (potential) bounds
        ss.low.list <- NULL
        for(m in 1:(M-numZero) )
        {
            ss.low.list[[m]] <- pump_sample_raw(
                d_m = d_m, MTP = MTP, typesample = typesample,
                MDES = MDES[m], J = J, K = K,
                target.power = target.power,
                nbar = nbar, Tbar = Tbar,
                alpha = alpha,
                two.tailed = two.tailed,
                numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, 
                numCovar.3 = numCovar.3,
                R2.1 = R2.1[m], R2.2 = R2.2[m], R2.3 = R2.3[m],
                ICC.2 = ICC.2[m], ICC.3 = ICC.3[m],
                omega.2 = omega.2[m], omega.3 = omega.3[m],
                warn.small = FALSE
            )
        }
        
        # Identify sample size for Bonferroni
        ss.high.list <- NULL
        for(m in 1:(M-numZero) )
        {
            ss.high.list[[m]] <- pump_sample_raw(
                d_m = d_m, MTP = MTP, typesample = typesample,
                MDES = MDES[m], J = J, K = K,
                target.power = target.power,
                nbar = nbar, Tbar = Tbar,
                alpha = alpha / M, # adjust alpha for BF
                two.tailed = two.tailed,
                numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, 
                numCovar.3 = numCovar.3,
                R2.1 = R2.1[m], R2.2 = R2.2[m], R2.3 = R2.3[m],
                ICC.2 = ICC.2[m], ICC.3 = ICC.3[m],
                omega.2 = omega.2[m], omega.3 = omega.3[m],
                warn.small = FALSE
            )
        }
        
    } else {
        # lower bound needs to be lower for min type power
        need_pow <- 1 - (1 - target.power)^(1/M)
        
        ss.low.list <- NULL
        for(m in 1:(M-numZero))
        {
            ss.low.list[[m]] <- pump_sample_raw(
                d_m = d_m, MTP = MTP, typesample = typesample,
                MDES = MDES[m], J = J, K = K,
                target.power = need_pow,
                nbar = nbar, Tbar = Tbar,
                alpha = alpha,
                two.tailed = two.tailed,
                numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, 
                numCovar.3 = numCovar.3,
                R2.1 = R2.1[m], R2.2 = R2.2[m], R2.3 = R2.3[m],
                ICC.2 = ICC.2[m], ICC.3 = ICC.3[m],
                omega.2 = omega.2[m], omega.3 = omega.3[m],
                warn.small = FALSE
            )
        }
        
        # higher bound needs to be higher for min type power (including complete)
        need_pow <- (target.power^(1/M))
        ss.high.list <- NULL
        for(m in 1:(M-numZero))
        {
            ss.high.list[[m]] <- pump_sample_raw(
                d_m = d_m, MTP = MTP, typesample = typesample,
                MDES = MDES[m], J = J, K = K,
                target.power = need_pow,
                nbar = nbar, Tbar = Tbar,
                alpha = alpha / M, # adjust alpha for BF
                two.tailed = two.tailed,
                numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, 
                numCovar.3 = numCovar.3,
                R2.1 = R2.1[m], R2.2 = R2.2[m], R2.3 = R2.3[m],
                ICC.2 = ICC.2[m], ICC.3 = ICC.3[m],
                omega.2 = omega.2[m], omega.3 = omega.3[m],
                warn.small = FALSE
            )
        }
    }
    
    ss.low.vals <- vapply(ss.low.list, function(x) x$ss, numeric(1))
    which.ss.low <- which.min(ss.low.vals)
    # check if everything is NA
    if(length(which.ss.low) > 0)
    {
        ss.low <- ss.low.list[[which.ss.low]]$ss
    } else
    {
        ss.low <- 1
    }
    
    ss.high.vals <- vapply(ss.high.list, function(x) x$ss, numeric(1))
    which.ss.high <- which.max(ss.high.vals)
    # check if everything is NA
    if(length(which.ss.high) > 0)
    {
        ss.high <- ss.high.list[[which.ss.high]]$ss
        default.max <- FALSE
    } else {
        if( typesample == 'nbar')
        {
            ss.high <- max_sample_size_nbar
        } else
        {
            ss.high <- max_sample_size_JK
        }
        default.max <- TRUE
    }
    
    # Done if Bonferroni is what we are looking for
    if (MTP == "BF" & pdef$indiv ) {
        ss.results <- data.frame(MTP, typesample, ss.high, target.power)
        colnames(ss.results) <- output.colnames
        
        if(default.max) {
            warning('Cannot achieve target power with given parameters.')
            ss.results <- data.frame(MTP, typesample, NA, target.power)
            colnames(ss.results) <- output.colnames
        }
        return( make.pumpresult( ss.results, tries = NULL,
                                 type = "sample", params.list = params.list,
                                 d_m = d_m,
                                 sample.level = typesample,
                                 exact = TRUE,
                                 power.params.list = pow_params) )
    }
    
    # Done if None is what we are looking for
    if (MTP == "None" ) {
        ss.results <- data.frame(MTP, typesample, ss.low, target.power)
        colnames(ss.results) <- output.colnames
        
        if(default.max)
        {
            warning('Cannot achieve target power with given parameters.')
            ss.results <- data.frame(MTP, typesample, NA, target.power)
            colnames(ss.results) <- output.colnames
        }
        
        return( make.pumpresult( ss.results, tries = NULL,
                                 type = "sample", params.list = params.list,
                                 d_m = d_m,
                                 sample.level = typesample,
                                 exact = TRUE,
                                 power.params.list = pow_params) )
    }
    
    if(default.max)
    {
        warning( "Using default max sample size for one end of initial bounds of search, so estimation may take more time.", call. = FALSE )
    }
    
    # search in the grid from min to max.
    test.pts <- optimize_power(
        d_m = d_m, search.type = typesample,
        MTP, target.power, power.definition, tol,
        start.low = ss.low, start.high = ss.high,
        MDES = MDES,
        J = J, K = K, nbar = nbar,
        M = M, numZero = numZero, Tbar = Tbar, 
        alpha = alpha, two.tailed = two.tailed,
        numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, 
        numCovar.3 = numCovar.3,
        R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, 
        ICC.2 = ICC.2, ICC.3 = ICC.3,
        rho = rho, omega.2 = omega.2, omega.3 = omega.3,
        B = B, parallel.WY.cores = parallel.WY.cores,
        max.steps = max.steps, 
        tnum = tnum, start.tnum = start.tnum, final.tnum = final.tnum,
        give.warnings = give.optimizer.warnings
    )
    
    # Assemble results
    # round up to get nice sufficient sample size.
    ss.results <- data.frame(
        MTP,
        typesample,
        ifelse(is.na(test.pts$pt[nrow(test.pts)]),
               NA,   # failed to find solution
               ceiling(test.pts$pt[nrow(test.pts)])),  
        test.pts$power[nrow(test.pts)]
    )
    colnames(ss.results) <- output.colnames
    
    
    # if it has converged, give notice about possible flatness
    if(is.finite(ss.results$`Sample.size`) && 
       test.pts$dx[[nrow(test.pts)]] < 0.005 ) {
        msg <- "Power curve is relatively flat. Other (smaller values) may have similar power.\nPlease refer to sample size vignette for interpretation."
        message(msg)
        flat <- TRUE
    } else {
        flat <- FALSE
    }
    
    return( make.pumpresult( ss.results, type = "sample", 
                             params.list = params.list,
                             d_m = d_m,
                             sample.level = typesample,
                             power.params.list = pow_params,
                             tries = test.pts,
                             flat = flat) )
}


