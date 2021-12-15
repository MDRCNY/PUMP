

#' List all the supported designs of the `pum` package.
#'
#' List all supported designs, with brief descriptions.
#'
#' @param comment = TRUE prints out description of each design or method.  FALSE does not.
#'
#' @export
pump_info <- function( comment = TRUE) {
    design <- tibble::tribble(
        ~ Code, ~PowerUp, ~ Comment, ~Params,
        # 1 level design
        "d1.1_m1c",     "n/a",
            "1 lvl, lvl 1 rand / constant impacts model",
            "R2.1",

        # 2 level designs, randomization at level 1
        "d2.1_m2fc",    "bira2_1c",
            "2 lvls, lvl 1 rand / fixed intercepts, constant impacts",
            "R2.1, ICC.2",
        "d2.1_m2ff",    "bira2_1f",
            "2 lvls, lvl 1 rand / fixed intercepts, fixed impacts",
            "R2.1, ICC.2",
        "d2.1_m2fr",    "bira2_1r",
            "2 lvls, lvl 1 rand / fixed intercepts, random impacts (FIRC)",
            "R2.1, ICC.2, omega.2",
        "d2.1_m2rr",    "n/a",
            "2 lvls, lvl 1 rand / random intercepts & impacts (RIRC)",
            "R2.1, ICC.2, omega.2",

        # 2 lvl design, rand at lvl 2
        "d2.2_m2rc",    "cra2_2r",
            "2 lvls, lvl 2 rand / random intercepts, constant impacts",
            "R2.1, R2.2, ICC.2",


        # 3 lvl design, rand at lvl 1
        "d3.1_m3rr2rr", "bira3_1r ",
            "3 lvls, lvl 1 rand / lvl 3 random intercepts, random impacts, lvl 2 random intercepts, random impacts",
            "R2.1, ICC.2, omega.2, ICC.3, omega.3",
        # 3 lvl design, rand at lvl 2
        "d3.2_m3ff2rc", "bcra3_2f",
            "3 lvls, lvl 2 rand / lvl 3 fixed intercepts, fixed impacts, lvl 2 random intercepts, constant impacts",
            "R2.1, R2.2, ICC.2, ICC.3",
        "d3.2_m3fc2rc", "n/a",
            "3 lvls, lvl 2 rand / lvl 3 fixed intercepts, constant impact, lvl 2 random intercepts, constant impact",
            "R2.1, R2.2, ICC.2, ICC.3",
        "d3.2_m3rr2rc", "bcra3_2r",
            "3 lvls, lvl 2 rand / lvl 3 random intercepts, random impacts, lvl 2 random intercepts, constant impacts",
            "R2.1, R2.2, ICC.2, ICC.3, omega.3",
        # 3 lvl design, rand at lvl 3
        "d3.3_m3rc2rc", "cra3_3r",
            "3 lvls, lvl 3 rand / lvl 3 random intercepts, constant impacts, lvl 2 random intercepts, constant impacts",
            "R2.1, R2.2, ICC.2, R2.3, ICC.3"
    )

    design <- tidyr::separate( design, .data$Code, into = c("Design", "Model"), remove = FALSE, sep = "_" )

    adjust <- tibble::tribble( ~ Method, ~ Comment,
                               "None", "No adjustment",
                               "Bonferroni", "The classic (and conservative) multiple testing correction",
                               "Holm", "Step down version of Bonferroni",
                               "BH", "Benjamini-Hochberg",
                               "WY-SS", "Westfall-Young, Single Step",
                               "WY-SD", "Westfall-Young, Step Down" )

    params <- tibble::tribble( ~ Parameter, ~ Description,
                               "nbar",       "the harmonic mean of the number of level 1 units per level 2 unit (students per school)",
                               "J",          "the number of level 2 units (schools)",
                               "K",          "the number of level 3 units (district)",
                               "Tbar",       "the proportion of units that are assigned to the treatment",
                               "numCovar.1", "number of Level 1 (individual) covariates",
                               "numCovar.2", "number of Level 2 (school) covariates",
                               "numCovar.3", "number of Level 3 (district) covariates",
                               "R2.1",       "percent of variation explained by Level 1 covariates",
                               "R2.2",       "percent of variation explained by Level 2 covariates",
                               "R2.3",       "percent of variation explained by Level 3 covariates",
                               "ICC.2",      "level 2 intraclass correlation",
                               "ICC.3",      "level 3 intraclass correlation",
                               "omega.2",    "ratio of variance of level 2 average impacts to variance of level 2 random intercepts",
                               "omega.3",    "ratio of variance of level 3 average impacts to variance of level 3 random intercepts"
    )

    if ( !comment ) {
        design$Comment <- NULL
        adjust$Comment <- NULL
    }

    list( Design = design, Adjustment = adjust, Parameters = params )
}




#' Return characteristics of a given design code
#'
#' See the pump_info method to get a list of supported designs.
#'
#' @param design String. Experimental design to parse.
#'
#' @return List of features including number of levels, level of randomization,
#'   etc.
#'
#' @family pump_info
#'
#' @examples
#' supported = pump_info()$Design
#' supported$Code[[4]]
#' parse_design( supported$Code[[4]] )
#'
#' @export
parse_design <- function( design ) {
    des <- stringr::str_split(design, "\\.|_")[[1]]
    nums <- readr::parse_number(des)
    levels <- nums[[1]]
    if ( levels == 3 ) {
        l3 <- substr( des[3], 3, 4)
        l3.p <- strsplit( l3, "" )[[1]]
        l2 <- substr( des[3], 6, 8 )
        l2.p <- strsplit( l2, "" )[[1]]
    } else if ( levels == 2 ) {
        l2 <- substr( des[3], 2, 4 )
        l2.p <- strsplit( l2, "" )[[1]]
        l3 <- NULL
        l3.p <- NULL
    } else {
        l2 <- NULL
        l2.p <- NULL
        l3 <- NULL
        l3.p <- NULL
    }

    FE.2 <- !is.na(l2) && substring( l2, 0, 1 ) == "f"
    FE.3 <- !is.na(l3) && substring( l3, 0, 1 ) == "f"

    list( levels = levels,
          rand_level = nums[[2]],
          model2 = l2,
          model2.p = l2.p,
          model3 = l3,
          model3.p = l3.p,
          FE.2 = FE.2,
          FE.3 = FE.3
    )
}






#' Computes Q_m, the standard error of the effect size estimate
#'
#' Function to calculate the theoretical true (unadjusted) standard error of the
#' ATE estimate for a given design and model, in effect size units.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param J scalar; the number of schools
#' @param K scalar; the number of districts
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that are assigned to the
#'   treatment
#' @param R2.1 scalar, or vector of length M; percent of variation explained by
#'   Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M; percent of variation explained by
#'   Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M; percent of variation explained by
#'   Level 3 covariates for each outcome
#' @param ICC.2 scalar, or vector of length M; school intraclass correlation
#' @param ICC.3 scalar, or vector of length M; district intraclass correlation
#' @param omega.2 scalar, or vector of length M; ratio of school effect size
#'   variability to random effects variability
#' @param omega.3 scalar, or vector of length M; ratio of district effect size
#'   variability to random effects variability
#'
#' @return Q_m, the standard error of the effect size estimate
#' @export

calc_SE <- function(design, J, K, nbar, Tbar, R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3) {

    if(design %in% c('d1.1_m1c'))
    {
        Q.m <- sqrt( ( (1 - R2.1) )  /(Tbar * (1-Tbar) * nbar) )
    } else if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
    {
        Q.m <- sqrt( ( (1 - ICC.2)*(1 - R2.1) ) / (Tbar * (1-Tbar) * J * nbar) )
    } else if (design == 'd2.1_m2fr' || design == 'd2.1_m2rr' )
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
    } else if (design == 'd3.2_m3ff2rc' || design == 'd3.2_m3fc2rc' )
    {
        Q.m <- sqrt( ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K) ) +
                    ( ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * K * nbar) ) )
    } else if (design == 'd3.2_m3rr2rc' )
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
#' @param validate whether or not to validate if output df is <= 0
#'
#' @return Degree of freedom for the design.
#'
#' @export

calc_df <- function(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3, validate = TRUE) {

    if(design == 'd1.1_m1c')
    {
        df <- J * nbar - numCovar.1 - 1
    } else if(design == 'd2.1_m2fc')
    {
        df <- J * (nbar - 1) - numCovar.1 - 1
    } else if (design == 'd2.1_m2ff')
    {
        df <- J * (nbar - 2) - numCovar.1
    } else if (design == 'd2.1_m2fr' || design == 'd2.1_m2rr' )
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
    } else if (design == 'd3.2_m3ff2rc' )
    {
        df <- K * (J - 2) - numCovar.2
    } else if (design == 'd3.2_m3fc2rc' ) {
        df <- K * (J - 1) - numCovar.2
    } else if (design == 'd3.2_m3rr2rc')
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





#' This function calculates needed nbar to achieve a given power (as represented by
#' a difference in t-statistic) for all implemented designs
#'
#' @inheritParams calc_J
#' @param J scalar; the number of schools
#'
#' @return nbar, the number of individuals needed, or NA if not possible given design
#' @export

calc_nbar <- function(design, MT = 2.8, MDES, J, K = NULL, Tbar, R2.1,
                      R2.2, ICC.2, omega.2,
                      R2.3 = NULL, ICC.3 = NULL, omega.3 = NULL ) {

    if(design %in% c('d1.1_m1c'))
    {
        numr <- (1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J
        nbar <- (MT/MDES)^2 * numr/denom
    } else if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
    {
        numr <- (1 - ICC.2) * (1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J
        nbar <- (MT/MDES)^2 * numr/denom
    } else if (design == 'd2.1_m2fr' || design == 'd2.1m2rr' )
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
        denom <- Tbar * (1 - Tbar) * J * K * ((MDES/MT)^2) - J * ICC.3 * (1 - R2.3)  - ICC.2 * (1 - R2.2)
        nbar <- numr / denom
    } else if (design == 'd3.2_m3ff2rc' || design == 'd3.2_m3fc2rc' )
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
#' @param MDES scalar; the MDES values for each outcome
#'
#' @return J, the number of schools needed
#' @export

calc_J <- function(
    design, MT = 2.8, MDES, K = NULL, nbar, Tbar,
    R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3
) {

    if(design %in% c('d1.1_m1c'))
    {
        numr <- (1 - R2.1)
        denom <- (Tbar * (1 - Tbar) * nbar)
        J <- (MT/MDES)^2 * numr/denom
    } else if(design %in% c('d2.1_m2fc', 'd2.1_m2ff'))
    {
        numr <- (1 - ICC.2) * (1 - R2.1)
        denom <- (Tbar * (1 - Tbar) * nbar)
        J <- (MT/MDES)^2 * numr/denom
    } else if (design == 'd2.1_m2fr' || design == 'd2.1m2rr' )
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
    } else if (design == 'd3.2_m3ff2rc' || design == 'd3.2_m3fc2rc' )
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
#' @param MDES scalar; the MDES value for all outcomes
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
calc_K <- function(design, MT, MDES, J, nbar, Tbar,
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
    } else if (design == 'd3.2_m3ff2rc' || design == 'd3.2_m3fc2rc' )
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



make_MDES_vector = function( MDES, M, numZero = NULL, verbose = TRUE ) {
    if( !is.null(numZero) ) {
        if( ( length(MDES) > 1 ) && ( numZero + length(MDES) != M ) )
        {
            stop('Please provide an MDES vector + numZero that add up to M.\n
             Example: MDES = c(0.1, 0.1), numZero = 3, M = 5.\n
             Assumed MDES vector = c(0.1, 0.1, 0, 0, 0)')
        }
        if( numZero >= M )
        {
            stop('numZero cannot be greater than or equal to M' )
        }
        if ( length(MDES) == 1 ) {
            MDES <- c(rep( MDES, M - numZero), rep(0, numZero) )
        } else {
            MDES <- c(MDES, rep(0, numZero))
        }
        if ( verbose ) {
            message('Assumed full MDES vector:', 'c(', paste(MDES, collapse = ', '), ')')
        }
    }

    if(length(MDES) != M)
    {
        if ( length(MDES) == 1 ) {
            MDES <- rep( MDES, M )
        } else {
            stop(paste('Please provide a vector of MDES values of length 1 or M. Current vector:',
                       MDES, 'M =', M))
        }
    }

    MDES
}


validate_MTP = function( MTP, power.call, M, multi.MTP.ok = FALSE ) {

    if( !multi.MTP.ok && length( MTP ) > 1 )
    {
        stop( 'Please provide only a single MTP procedure.' )
    }

    if( !is.null( MTP ) && any( MTP == "raw" ) ) {
        MTP[ MTP == "raw" ] = "None"
    }

    if(M == 1)
    {
        if ( !is.null( MTP ) && (MTP != "None" ) )
        {
            warning("Multiple testing corrections are not needed when M = 1.")
        }
        MTP <- "None"
    } else if( power.call && (is.null(MTP) || MTP == 'None')) {
        stop('Please provide a multiple test procedure (MTP).')
    }

    chk = MTP %in% pump_info()$Adjustment$Method

    if( ! all( chk ) ) {
        if ( length( MTP ) > 1 ) {
            msg = sprintf( 'You have at least one invalid MTP: %s',
                       paste( "'", MTP[!chk], "'", sep = "", collapse = ", " ) )
        } else {
            msg = sprintf( '"%s" is an invalid MTP.', MTP )
        }
        stop( msg )
    }

    return( MTP )
}



#' Validates user inputs
#'
#' This functions takes in a list of user inputs. Depending on the inputs,
#' it produces errors or warnings, and at times modifies inputs if necessary.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param params.list a list of parameters input by a user
#' @param power.call flag for power estimation
#' @param ss.call flag for sample size estimation
#' @param mdes.call flag for MDES estimation
#'
#' @return params.list
#'
validate_inputs <- function( design, params.list,
                             power.call = FALSE,
                             mdes.call = FALSE,
                             ss.call = FALSE,
                             verbose = TRUE,
                             multi.MTP.ok = FALSE )
{

    #-------------------------------------------------------#
    # basic checks of inputs
    #-------------------------------------------------------#

    # allow either supported design names or PowerUp equivalents
    info <- pump_info()
    if(!(design %in% info$Design$Code))
    {
        if(design %in% info$Design$PowerUp)
        {
            design <- info$Design$Code[info$Design$PowerUp == design]
        } else
        {
            stop('Invalid design.')
        }
    }

    par_design <- parse_design(design)

    params.list$MTP = validate_MTP( MTP = params.list$MTP,
                                    power.call= power.call,
                                    M = params.list$M,
                                    multi.MTP.ok = multi.MTP.ok )

    #-------------------------------------------------------#
    # Drop inputs that we can ignore, depending on model
    #-------------------------------------------------------#

    # check for three-level parameters when working with two level models
    if( par_design$levels <= 2 )
    {

        if( ( !is.null(params.list$K) && params.list$K > 1 ) |
            ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
            ( !is.null(params.list$R2.3)) && any( params.list$R2.3 > 0 ) |
            ( !is.null(params.list$omega.3) && any(params.list$omega.3 > 0 ) ) )
        {
            warning('The following parameters are only valid for three-level designs, and will be ignored:\n
              K, numCovar.3, R2.3, ICC.3, omega.3')
            params.list$K <- NULL
            params.list$numCovar.3 <- NULL
            params.list$R2.3 <- NULL
            params.list$ICC.3 <- NULL
            params.list$omega.3 <- NULL
        }
    }

    if( par_design$levels == 3 && par_design$FE.3 )
    {
        if( ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
            ( !is.null(params.list$R2.3) && any( params.list$R2.3 > 0 ) ) )
        {
            warning('The following parameters are not valid for fixed effect designs, and will be ignored:\n
              numCovar.3, R2.3')
            params.list$numCovar.3 <- NULL
            params.list$R2.3 <- NULL
        }
    }



    #-------------------------------------------------------#
    # MDES
    #-------------------------------------------------------#

    if ( mdes.call ) {
        if ( !is.null( params.list$MDES ) ) {
            stop( "You cannot provide MDES to pump_mdes()" )
        }
        if ( !is.null( params.list$numZero ) && params.list$numZero >= params.list$M ) {
            stop( sprintf( "You cannot specify %s zeros with %s outcomes", params.list$numZero, params.list$M ) )
        }
    } else {
        params.list$MDES = make_MDES_vector( params.list$MDES, params.list$M, params.list$numZero,
                                             verbose = verbose )
    }


    #---------------------------------------------------------------#
    # convert all params from scalar to vector, if they are non-null
    #---------------------------------------------------------------#
    if(!(length(params.list$R2.1) %in% c(1, params.list$M)))
    {
        stop("R2.1: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$R2.1) == 1)
    {
        params.list$R2.1 <- rep(params.list$R2.1, params.list$M)
    }

    if(!(length(params.list$R2.2) %in% c(0, 1, params.list$M)))
    {
        stop("R2.2: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$R2.2) == 1)
    {
        params.list$R2.2 <- rep(params.list$R2.2, params.list$M)
    }

    if(!(length(params.list$R2.3) %in% c(0, 1, params.list$M)))
    {
        stop("R2.3: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$R2.3) == 1)
    {
        params.list$R2.3 <- rep(params.list$R2.3, params.list$M)
    }

    if(!(length(params.list$ICC.2) %in% c(0, 1, params.list$M)))
    {
        stop("ICC.2: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$ICC.2) == 1)
    {
        params.list$ICC.2 <- rep(params.list$ICC.2, params.list$M)
    }

    if(!(length(params.list$ICC.3) %in% c(0, 1, params.list$M)))
    {
        stop("ICC.3: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$ICC.3) == 1)
    {
        params.list$ICC.3 <- rep(params.list$ICC.3, params.list$M)
    }

    if(!(length(params.list$omega.2) %in% c(0, 1, params.list$M)))
    {
        stop("omega.2: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$omega.2) == 1)
    {
        params.list$omega.2 <- rep(params.list$omega.2, params.list$M)
    }

    if(!(length(params.list$omega.3) %in% c(0, 1, params.list$M)))
    {
        stop("omega.3: Please provide a scalar parameter or a vector of length M.")
    }
    if(length(params.list$omega.3) == 1)
    {
        params.list$omega.3 <- rep(params.list$omega.3, params.list$M)
    }




    #-------------------------------------------------------#
    # Basic checks of data parameters
    #-------------------------------------------------------#

    if( (!is.null( params.list$K ) && params.list$K <= 0) |
        ( !is.null( params.list$J) && params.list$J <= 0) |
        params.list$nbar <= 0)
    {
        stop('Provided values of J, K, and/or nbar need to be positive.')
    }

    if( params.list$numCovar.1 < 0 |
        ( !is.null( params.list$numCovar.2 )  && params.list$numCovar.2 < 0  ) |
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
        stop('Please provide alpha as a probability between 0 and 1')
    }

    if(any(params.list$R2.1 > 1) | any(params.list$R2.1 < 0) |
       any(params.list$R2.2 > 1) | any(params.list$R2.2 < 0) |
       any(params.list$R2.3 > 1) | any(params.list$R2.3 < 0))
    {
        stop('Please provide R2 as a probability between 0 and 1')
    }

    if(any(params.list$omega.2 < 0) | (!is.null(params.list$omega.3) && any(params.list$omega.3 < 0)))
    {
        stop('Please provide a non-negative value of Omega')
    }


    #-------------------------------------------------------#
    # check for inconsistent user inputs
    #-------------------------------------------------------#

    # one level models
    if( par_design$levels == 1 ) {
        if ( !is.null( params.list$J ) && params.list$J != 1 ) {
            stop( "Can't have multiple units at 2nd level for the d1.1_m2cc design" )
        }

        #  if ( !is.null( params.list$numCovar.2 ) || !is.null( params.list$numCovar.3 ) || !is.null( params.list$R2.2 ) || !is.null( params.list$ICC.2 ) ) {
        #   stop( "Can't have numCovar.2, numCovar.3, R2.2, or ICC.2 for d1.1_m2cc design (single level design" )
        #  }

    }

    # two level models
    if( par_design$levels <= 2 )
    {
        if ( par_design$levels == 2 & params.list$J == 1 )
        {
            warning('Two level design with single unit at level 2')
        }
    }


    # three level models
    if( par_design$levels == 3 )
    {
        if(is.null(params.list$K) || params.list$K < 1 )
        {
            stop('You must specify K, with K >= 1 (number of units at level 3) for three-level designs' )
        }
        if(params.list$K == 1)
        {
            warning('Running a 3-level model with K = 1')
        }
    }

    # three level models, continued.
    if( par_design$levels == 3 && !par_design$FE.3 )
    {
        if( is.null(params.list$numCovar.3) | is.null(params.list$R2.3))
        {
            stop('You must specify both numCovar.3 and R2.3 for three-level designs with random effects')
        }
    }

    # ICC
    if(!is.null(params.list$ICC.2) && !is.null(params.list$ICC.3) && any(params.list$ICC.2 + params.list$ICC.3 > 1))
    {
        stop('ICC.2 + ICC.3 must be <= 1')
    }

    # constant treatment effects models: level 2
    if( par_design$levels >= 2 && par_design$model2.p[[2]] == 'c' )
    {
        if(any(params.list$omega.2 > 0))
        {
            warning('Omega is assumed to be 0 for constant treatment effects models. Ignoring input omega.2 value')
            params.list$omega.2 <- 0
        }
    }

    # constant treatment effects models: level 3
    if( par_design$levels == 3 && par_design$model3 %in% c('rc', 'cc'))
    {
        if(!is.null(params.list$omega.3) && any(params.list$omega.3 > 0))
        {
            warning('Omega is assumed to be 0 for constant treatment effects models. Ignoring input omega.3 value')
            params.list$omega.3 <- 0
        }
    }

    # If 3 level model allows treatment variation, we need omega
    if( par_design$levels == 3 && par_design$model3.p[2] != 'c' )
    {
        if(is.null(params.list$omega.3))
        {
            stop('Omega.3 is required for this design.')
        }
    }

    # NOTE: Used to be this:     design %in% c('d3.1_m3rr2rr', 'd3.2_m3rr2rc'))
    # But shouldn't all three level designs need ICC.3 (possibly = 0, of course).
    if( par_design$levels == 3 )
    {
        if(is.null(params.list$ICC.3))
        {
            stop('ICC.3 is required for this design.')
        }
    }

    # number covariates
    if(!is.null( params.list$R2.1) && any(params.list$R2.1 != 0) && params.list$numCovar.1 == 0)
    {
        warning('If nonzero R2 (R2.1, level 1), at least one covariate is assumed. Setting numCovar.1 = 1')
        params.list$numCovar.1 <- 1
    }
    if(!is.null( params.list$R2.2) && any(params.list$R2.2 != 0) && params.list$numCovar.2 == 0)
    {
        warning('If nonzero R2 (R2.2, level 2), at least one covariate is assumed. Setting numCovar.2 = 1')
        params.list$numCovar.2 <- 1
    }
    if(!is.null( params.list$R2.3) && any(params.list$R2.3 != 0) && params.list$numCovar.3 == 0)
    {
        warning('If nonzero R2 (R2.3, level 3), at least one covariate is assumed. Setting numCovar.3 = 1')
        params.list$numCovar.3 <- 1
    }

    #-------------------------------------------------------#
    #  rho
    #-------------------------------------------------------#
    if(is.null(params.list$rho.matrix) & is.null(params.list$rho))
    {
        stop( sprintf( 'Please provide either a %d x %d rho.matrix or default scalar rho.',
                       params.list$M, params.list$M ) )
    }

    if(!is.null(params.list$rho.matrix))
    {
        if(nrow(params.list$rho.matrix) != params.list$M | ncol(params.list$rho.matrix) != params.list$M)
        {
            stop('Correlation matrix of invalid dimensions. Please provide valid correlation matrix.')
        }
    } else {
        if(params.list$rho > 1 | params.list$rho < -1)
        {
            stop('Please provide rho as a correlation between -1 and 1')
        }
    }

    return(params.list)

}



# Stolen from development purrr
silently <- function(.f, otherwise = NULL) {
    .f <- as_mapper(.f)
    function(...) {
        ret <-
            purrr:::capture_output(
                purrr:::capture_error(.f(...), otherwise, quiet=TRUE)
            )
        # reformat output to an un-nested list
        list(
            result = ret$result$result,
            output = ret$output,
            messages = ret$messages,
            warnings = ret$warnings,
            error = ret$result$error
        )
    }
}


# Not yet implemented
# #' @param grid Flag of whether call is a grid call or non-grid call.


#' Check user inputs
#'
#' This functions takes in a list of user inputs and checks them for validity,
#' producing a mix of errors or warnings.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param call String denoting intended pump_power, pump_mdes, or pump_sample
#'   call.
#' @param verbose Keep and return smaller messages or not
#' @param multi.MTP.ok TRUE/FALSE on whether multiple MTP is allowed for check.
#' @param ... The arguments to be passed to given call.
#'
#' @return list of two things: "ok", a TRUE/FALSE of whether this is a valid
#'   call.  "messages", a list of all messages, warnings, and notes that would
#'   be generated by the input parameter choices.  In the case of an error, the
#'   messages will have the error as the first string in the list.
#' @import purrr
#' @export
check_pump_call = function( design,
                            call = c( "power", "mdes", "sample" ),
                            #grid = FALSE,
                            verbose = TRUE,
                            multi.MTP.ok = TRUE,
                            ... ) {

    params = list( ... )
    call = match.arg(call)
    power.call = call == "power"
    mdes.call = call == "mdes"
    ss.call = call == "sample"
    params$design = design


    quiet_validate_inputs = silently( validate_inputs )
    cres = quiet_validate_inputs( design = design, params.list = params,
                                 power.call = power.call, mdes.call = mdes.call, ss.call = ss.call,
                                 verbose = verbose,
                                 multi.MTP.ok = multi.MTP.ok )

    res = list()
    res$ok = is.null( cres$error )

    messages = c()
    if ( !res$ok ) {
        messages = cres$error$message
    }
    messages = c( messages, cres$warnings, cres$messages )
    if ( cres$output != "" ) {
        messages = c( messages, cres$output )
    }
    res$messages = messages
    return( res )
}


# print out results cleanly
scat <- function( str, ... ) {
    cat( sprintf( str, ... ) )
}


