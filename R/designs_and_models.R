

#' List all the supported designs of the `pum` package.
#'
#' List all supported designs, with brief descriptions.
#'
#' @param comment = TRUE prints out description of each design or method.  FALSE does not.
#'
#' @export
pump_info <- function( comment = TRUE) {
    design <- tibble::tribble(
        ~ Code, ~PowerUp, ~ Comment,
        # 1 level design
        "d1.1_m1c",     "n/a",            "1 lvl, lvl 1 rand / constant impacts model",
        # 2 level designs, randomization at level 1
        "d2.1_m2fc",    "blocked_i1_2c",  "2 lvls, lvl 1 rand / fixed intercepts, constant impacts",
        "d2.1_m2ff",    "blocked_i1_2f",  "2 lvls, lvl 1 rand / fixed intercepts, fixed impacts",
        "d2.1_m2fr",    "blocked_i1_2r",  "2 lvls, lvl 1 rand / fixed intercepts, random impacts",
        # 2 lvl design, rand at lvl 2
        "d2.2_m2rc",    "simple_2c_2r",   "2 lvls, lvl 2 rand / random intercepts, constant impacts",
        # 3 lvl design, rand at lvl 1
        "d3.1_m3rr2rr", "blocked_i1_3r",  "3 lvls, lvl 1 rand / lvl 3 random intercepts, random impacts, lvl 2 random intercepts, random impacts",
        # 3 lvl design, rand at lvl 2
        "d3.2_m3ff2rc", "blocked_c2_3f",  "3 lvls, lvl 2 rand / lvl 3 fixed intercepts, fixed impacts, lvl 2 random intercepts, constant impacts",
        "d3.2_m3fc2rc", "n/a",            "3 lvls, lvl 2 rand / lvl 3 fixed intercepts, constant impact, lvl 2 random intercepts, constant impact",
        "d3.2_m3rr2rc", "blocked_c2_3r",  "3 lvls, lvl 2 rand / lvl 3 random intercepts, random impacts, lvl 2 random intercepts, constant impacts",
        # 3 lvl design, rand at lvl 3
        "d3.3_m3rc2rc", "simple_c3_3r",   "3 lvls, lvl 3 rand / lvl 3 random intercepts, constant impacts, lvl 2 random intercepts, constant impacts"
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

calc.Q.m <- function(design, J, K, nbar, Tbar, R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3) {

    if(design %in% c('d1.1_m1c'))
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

calc.df <- function(design, J, K, nbar, numCovar.1, numCovar.2, numCovar.3, validate = TRUE) {

    if(design == 'd1.1_m1c')
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
#' @inheritParams calc.J
#' @param J scalar; the number of schools
#'
#' @return nbar, the number of individuals needed, or NA if not possible given design
#' @export

calc.nbar <- function(design, MT = 2.8, MDES, J, K = NULL, Tbar, R2.1,
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

calc.J <- function(
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



