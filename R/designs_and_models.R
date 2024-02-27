#' @title Provides details about supported package features (core function)
#'
#' @description List user options:
#' designs and models (d_m), including what
#' parameters are relevant for each context;
#' multiple testing procedures;
#' types of power;
#' design and model parameters.
#'
#' @seealso For more detailed information about user choices,
#' see the manuscript <doi:10.18637/jss.v108.i06>,
#' which includes a detailed Technical Appendix
#' including information about the designs and models
#' and parameters.
#'
#' @param topic string; what kind of info. One of:
#' all, context, adjustment, power, parameters.
#' @param comment logical; prints out
#' long description of
#' each design and method.
#'
#' @return list; a list of data frames with
#' information about each topic.
#'
#' @export
pump_info <- function(
    topic = c( "all", "context", "adjustment", "power", "parameters" ),
    comment = TRUE
) {
  topic <- match.arg( topic )

    context <- tibble::tribble(
        ~d_m, ~PowerUp, ~Params, ~ Comment,
        # 1 level design
        "d1.1_m1c",     "n/a",
            "R2.1",
            "1 lvl, lvl 1 rand / constant impacts model",

        # 2 level designs, randomization at level 1
        "d2.1_m2fc",    "bira2_1c",
            "R2.1, ICC.2",
            "2 lvls, lvl 1 rand / fixed intercepts, constant impacts",

        "d2.1_m2ff",    "bira2_1f",
            "R2.1, ICC.2",
            "2 lvls, lvl 1 rand / fixed intercepts, fixed impacts",

        "d2.1_m2fr",    "bira2_1r",
            "R2.1, ICC.2, omega.2",
            "2 lvls, lvl 1 rand / fixed intercepts, random impacts (FIRC)",

        "d2.1_m2rr",    "n/a",
            "R2.1, ICC.2, omega.2",
            "2 lvls, lvl 1 rand / random intercepts & impacts (RIRC)",

        # 2 lvl design, rand at lvl 2
        "d2.2_m2rc",    "cra2_2r",
            "R2.1, R2.2, ICC.2",
            "2 lvls, lvl 2 rand / random intercepts, constant impacts",

        # 3 lvl design, rand at lvl 1
        "d3.1_m3rr2rr", "bira3_1r ",
            "R2.1, ICC.2, omega.2, ICC.3, omega.3",
            "3 lvls, lvl 1 rand /
             lvl 3 random intercepts, random impacts,
             lvl 2 random intercepts, random impacts",

        "d3.1_m3ff2rr", "n/a",
        "R2.1, ICC.2, omega.2, ICC.3",
        "3 lvls, lvl 1 rand /
             lvl 3 fixed intercepts, fixed impacts,
             lvl 2 random intercepts, random impacts",

        # 3 lvl design, rand at lvl 2
        "d3.2_m3ff2rc", "bcra3_2f",
           "R2.1, R2.2, ICC.2, ICC.3",
            "3 lvls, lvl 2 rand /
            lvl 3 fixed intercepts, fixed impacts,
            lvl 2 random intercepts, constant impacts",

        "d3.2_m3fc2rc", "n/a",
            "R2.1, R2.2, ICC.2, ICC.3",
            "3 lvls, lvl 2 rand /
             lvl 3 fixed intercepts, constant impact,
             lvl 2 random intercepts, constant impact",

        "d3.2_m3rr2rc", "bcra3_2r",
            "R2.1, R2.2, ICC.2, ICC.3, omega.3",
            "3 lvls, lvl 2 rand /
             lvl 3 random intercepts, random impacts,
             lvl 2 random intercepts, constant impacts",

        # 3 lvl design, rand at lvl 3
        "d3.3_m3rc2rc", "cra3_3r",
            "R2.1, R2.2, ICC.2, R2.3, ICC.3",
            "3 lvls, lvl 3 rand /
             lvl 3 random intercepts, constant impacts,
             lvl 2 random intercepts, constant impacts",

    )


    context <- tidyr::separate(context, .data$d_m,
                               into = c("Design", "Model"),
                               remove = FALSE, sep = "_" )

    power <- tibble::tribble(
        ~ Definition, ~ Comment,
        "D*indiv",    "individual power for each nonzero outcome:
                       D1indiv, D2indiv, etc.",
        "indiv.mean", "mean across individual powers",
        "min*",       "probability to detect at least * outcomes:
                       min1, min2, etc. Up to number of
                       nonzero outcomes",
        "complete",   "probability to detect all outcomes;
                       NA if any outcomes are assumed to be zero"
    )


    adjust <- tibble::tribble( ~ Method, ~ Comment,
                               "None",  "No adjustment",
                               "BF",    "Bonferroni",
                               "HO",    "Holm",
                               "BH",    "Benjamini-Hochberg",
                               "WY-SS", "Westfall-Young, Single Step",
                               "WY-SD", "Westfall-Young, Step Down" )

    params <- tibble::tribble( ~ Parameter,  ~ Description,
      "nbar",       "scalar; harmonic mean of number of level 1 units per
                                level 2 unit (students per school)",
      "J",          "scalar; harmonic mean of number of level 2
                      units per level 3 unit (schools per district)",
      "K",          "scalar; number of level 3 units (districts)",
      "Tbar",       "scalar; proportion of units assigned to treatment",
      "numCovar.1", "scalar; number of level 1 (individual) covariates",
      "numCovar.2", "scalar; number of level 2 (school) covariates",
      "numCovar.3", "scalar; number of level 3 (district) covariates",
      "R2.1",       "scalar/vector; percent of level 1 variation explained by
                    covariates",
      "R2.2",       "scalar/vector; percent of level 2 variation explained by
                    covariates",
      "R2.3",       "scalar/vector; percent of level 3 variation explained by
                    covariates",
      "ICC.2",      "scalar/vector; level 2 intraclass correlation",
      "ICC.3",      "scalar/vector; level 3 intraclass correlation",
      "omega.2",    "scalar/vector; ratio of variance of level 2 average
                    impacts to level 2 random intercepts",
      "omega.3",    "scalar/vector; ratio of variance of level 3 average
                    impacts to level 3 random intercepts",
      "rho",        "scalar; correlation between all pairs of
                    test statistics",
      "rho.matrix", "matrix; full matrix of correlations between
                    test statistics"
    )

    if ( !comment ) {
        context$Comment <- NULL
        adjust$Comment <- NULL
    }


    res <- list( Context = context,
          Adjustment = adjust,
          Power = power,
          Parameters = params)
    if ( topic != "all" ) {
      names(res) <- tolower(names(res))
      return( res[[topic]] )
    } else {
      return( res)
    }

}




parse_design <- function(d_m) {
    pattern <- "^d([0-9])\\.([0-9])"   # pattern to match

    if (grepl(pattern, d_m)) {
        match <- regmatches(d_m, regexec(pattern, d_m))
        num1 <- match[[1]][[2]]
        num2 <- match[[1]][[3]]
        return( list( levels = num1,
                      rand_level = num2,
                      design = match[[1]][[1]]) )
    }
    return( NULL )

}


#' @title Return characteristics of a given context/d_m code (support function)
#'
#' @description Returns number of levels and model at each level.
#' See pump_info()$Context to get a list of supported d_ms.
#'
#' @param d_m string; context to parse.
#'
#' @return list; list of features including number of levels,
#' level of randomization, etc.
#'
#' @family pump_info
#'
#' @examples
#' supported <- pump_info(comment = FALSE)$Context
#' parse_d_m( supported$d_m[4] )
#'
#' @export
parse_d_m <- function(d_m) {
    des <- stringr::str_split(d_m, "\\.|_")[[1]]
    nums <- readr::parse_number(des)
    levels <- nums[[1]]
    if ( levels == 3 ) {
        l3 <- substr( des[3], 3, 4)
        l3.p <- strsplit( l3, "" )[[1]]
        l2 <- substr( des[3], 6, 7 )
        l2.p <- strsplit( l2, "" )[[1]]
    } else if ( levels == 2 ) {
        l2 <- substr( des[3], 3, 4 )
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

    return( list( levels = levels,
          rand_level = nums[[2]],
          model2 = l2,
          model2.p = l2.p,
          model3 = l3,
          model3.p = l3.p,
          FE.2 = FE.2,
          FE.3 = FE.3,
          design = paste0( "d", levels, ".", nums[[2]] )
    ) )
}






#' Computes Q_m, the standard error of the effect size estimate
#'
#' Function to calculate the theoretical true (unadjusted) standard
#' error of the ATE estimate for a given d_m and model, in effect size units.
#'
#' @param d_m a single RCT d_m (see list/naming convention).
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
#' @return vector; the standard error of the effect size estimate
#' @keywords internal
calc_SE <- function(d_m, J, K, nbar, Tbar,
                    R2.1, R2.2, R2.3, ICC.2, ICC.3,
                    omega.2, omega.3) {
  if (d_m == "d1.1_m1c") {
    Q.m <- sqrt(
      ((1 - R2.1)) / (Tbar * (1 - Tbar) * nbar)
    )
  } else if (d_m %in% c("d2.1_m2fc", "d2.1_m2ff")) {
    Q.m <- sqrt(
      ((1 - ICC.2) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar)
    )
  } else if (d_m %in% c("d2.1_m2fr", "d2.1_m2rr")) {
    Q.m <- sqrt(
      (ICC.2 * omega.2) / J +
        ((1 - ICC.2) * (1 - R2.1)) /
          (Tbar * (1 - Tbar) * J * nbar)
    )
  } else if (d_m == "d3.1_m3rr2rr") {
    Q.m <- sqrt(
      (ICC.3 * omega.3) / K +
        (ICC.2 * omega.2) / (J * K) +
        ((1 - ICC.2 - ICC.3) * (1 - R2.1)) /
          (Tbar * (1 - Tbar) * J * K * nbar)
    )
  } else if (d_m == "d3.1_m3ff2rr") {
    Q.m <- sqrt(
      (ICC.2 * omega.2) / (J * K) +
        ((1 - ICC.2 - ICC.3) * (1 - R2.1)) /
          (Tbar * (1 - Tbar) * J * K * nbar)
    )
  } else if (d_m == "d2.2_m2rc") {
    Q.m <- sqrt(
      (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
        (1 - ICC.2) * (1 - R2.1) /
          (Tbar * (1 - Tbar) * J * nbar)
    )
  } else if (d_m == "d3.3_m3rc2rc") {
    Q.m <- sqrt(
      (ICC.3 * (1 - R2.3)) / (Tbar * (1 - Tbar) * K) +
        (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K) +
        ((1 - ICC.2 - ICC.3) * (1 - R2.1)) /
          (Tbar * (1 - Tbar) * J * K * nbar)
    )
  } else if (d_m == "d3.2_m3ff2rc" || d_m == "d3.2_m3fc2rc") {
    Q.m <- sqrt(
      ((ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K)) +
        (((1 - ICC.2 - ICC.3) * (1 - R2.1)) /
          (Tbar * (1 - Tbar) * J * K * nbar))
    )
  } else if (d_m == "d3.2_m3rr2rc") {
    Q.m <- sqrt(
      ((ICC.3 * omega.3) / K) +
        ((ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J * K)) +
        (((1 - ICC.2 - ICC.3) * (1 - R2.1)) /
          (Tbar * (1 - Tbar) * J * K * nbar))
    )
  } else {
    stop(paste("d_m not implemented:", d_m))
  }
  return(Q.m)
}


#' @title Calculate degrees of freedom (support function)
#'
#' @description Given sample sizes, return the used degrees of freedom
#' (frequently conservative) for the design and model.
#'
#' @inheritParams pump_power
#' @param validate logical; whether or not to validate
#' if output df is <= 0.
#'
#' @return scalar; degrees of freedom for the context.
#' @export
calc_df <- function(d_m, J, K, nbar,
                    numCovar.1, numCovar.2, numCovar.3,
                    validate = TRUE) {

    if (d_m == 'd1.1_m1c') {
        df <- nbar - numCovar.1 - 2
    } else if (d_m == 'd2.1_m2fc') {
        df <- J * (nbar - 1) - numCovar.1 - 1
    } else if (d_m == 'd2.1_m2ff') {
        df <- J * (nbar - 2) - numCovar.1
    } else if (d_m == 'd2.1_m2fr' || d_m == 'd2.1_m2rr' ) {
        df <- J - numCovar.1 - 1
    } else if (d_m == 'd3.1_m3rr2rr')
    {
        df <- K - 1
    } else if (d_m == 'd3.1_m3ff2rr')
    {
        df <- K * (J - 1) - 1
    } else if (d_m == 'd2.2_m2rc')
    {
        df <- J - numCovar.1 - 2
    } else if (d_m == 'd3.3_m3rc2rc')
    {
        df <- K - numCovar.3 - 2
    } else if (d_m == 'd3.2_m3ff2rc' )
    {
        df <- K * (J - 2) - numCovar.2
    } else if (d_m == 'd3.2_m3fc2rc' )
    {
        df <- K * (J - 1) - numCovar.2
    } else if (d_m == 'd3.2_m3rr2rc')
    {
        df <- K - 1
    } else
    {
        stop(paste('d_m not implemented:', d_m))
    }

    if (validate & df <= 0)
    {
        stop('Invalid d_m parameters resulting in nonpositive degrees of freedom')
    }

    return(df)
}





#' This function calculates needed nbar to achieve a given power
#'
#' @inheritParams calc_J
#' @param J scalar; the number of schools
#'
#' @return nbar, the number of individuals needed,
#' or NA if not possible given design
#' @keywords internal
calc_nbar <- function(d_m, MT = 2.8, MDES,
                      J = NULL, K = NULL, Tbar, R2.1,
                      R2.2 = NULL, ICC.2 = NULL, omega.2 = NULL,
                      R2.3 = NULL, ICC.3 = NULL, omega.3 = NULL
) {

    if (d_m %in% c('d1.1_m1c'))
    {
        numr <- (1 - R2.1)
        denom <- Tbar * (1 - Tbar)
        nbar <- (MT/MDES)^2 * numr/denom
    } else if (d_m %in% c('d2.1_m2fc', 'd2.1_m2ff'))
    {
        numr <- (1 - ICC.2) * (1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J
        nbar <- (MT/MDES)^2 * numr/denom
    } else if (d_m == 'd2.1_m2fr' || d_m == 'd2.1m2rr' )
    {
        numr <- (1 - ICC.2) * (1 - R2.1)
        denom <- J * ((MDES / MT)^2) - ICC.2 * omega.2
        nbar <- numr / (Tbar * (1 - Tbar) * denom)
    } else if (d_m == 'd3.1_m3rr2rr') {
        numr <- (1 - ICC.2 - ICC.3) * (1 - R2.1)
        denom <- J * K * ((MDES / MT)^2) -
            J * ICC.3 * omega.3 - ICC.2 * omega.2
        nbar <- numr / ( Tbar * (1 - Tbar) * denom )
    } else if (d_m == 'd3.1_m3ff2rr') {
        numr <- (1 - ICC.2 - ICC.3) * (1 - R2.1)
        denom <- J * K * ((MDES / MT)^2) - ICC.2 * omega.2
        nbar <- numr / ( Tbar * (1 - Tbar) * denom )
    } else if (d_m == 'd2.2_m2rc')
    {
        numr <- (1 - ICC.2)*(1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J * ((MDES/MT)^2) - ICC.2 * (1 - R2.2)
        nbar <- numr / denom
    } else if (d_m == 'd3.3_m3rc2rc')
    {
        numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J * K *
          ((MDES/MT)^2) - J * ICC.3 * (1 - R2.3)  - ICC.2 * (1 - R2.2)
        nbar <- numr / denom
    } else if (d_m == 'd3.2_m3ff2rc' || d_m == 'd3.2_m3fc2rc' )
    {
        numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J * K *
          ((MDES/MT)^2) - ICC.2 * (1 - R2.2)
        nbar <- numr / denom
    } else if (d_m == 'd3.2_m3rr2rc')
    {
        numr <- (1 - ICC.2 - ICC.3)*(1 - R2.1)
        denom <- Tbar * (1 - Tbar) * J *
          ( K * ((MDES/MT)^2) - ICC.3 * omega.3 ) - ICC.2 * (1 - R2.2)
        nbar <- numr / denom
    } else
    {
        stop(paste('d_m not implemented:', d_m))
    }
    nbar <- ifelse( is.na( nbar ) || nbar < 0, NA, nbar )
    return( nbar )
}



#' This function calculates needed J to achieve a given (unadjusted) power
#'
#' @inheritParams pump_power
#'
#' @param d_m a single RCT design (see list/naming convention)
#' @param MT Number of approximate effect-size unit SEs
#' (adjusted for degrees of freedom issues) that the MDES
#' needs to be to achieve desired power.  E.g., 2.8 for normal theory.
#' @param MDES scalar; the MDES values for each outcome
#'
#' @return J, the number of schools needed
#' @keywords internal
calc_J <- function(
    d_m, MT = 2.8, MDES, K = NULL, nbar, Tbar,
    R2.1, R2.2, R2.3, ICC.2, ICC.3, omega.2, omega.3
) {

    if (d_m %in% c('d1.1_m1c'))
    {
        numr <- (1 - R2.1)
        denom <- (Tbar * (1 - Tbar) * nbar)
        J <- (MT/MDES)^2 * numr/denom
    } else if (d_m %in% c('d2.1_m2fc', 'd2.1_m2ff'))
    {
        numr <- (1 - ICC.2) * (1 - R2.1)
        denom <- (Tbar * (1 - Tbar) * nbar)
        J <- (MT/MDES)^2 * numr/denom
    } else if (d_m == 'd2.1_m2fr' || d_m == 'd2.1m2rr' )
    {
        numr <- (1 - ICC.2) * (1 - R2.1)
        denom <- (Tbar * (1 - Tbar) * nbar)
        J <- (MT/MDES)^2 * ( (ICC.2 * omega.2) + numr / denom)
    } else if (d_m == 'd3.1_m3rr2rr')
    {
        numr <- (1 - ICC.2 - ICC.3 ) * (1 - R2.1) +
            Tbar * (1 - Tbar) * nbar * ICC.2 * omega.2
        denom <- K * (MDES/MT)^2 - ICC.3 * omega.3
        J <- (1 / (Tbar * (1 - Tbar) * nbar)) * numr/denom
    } else if (d_m == 'd3.1_m3ff2rr')
    {
        Q = (MT/MDES)^2
        tm1 = (ICC.2 * omega.2) / K
        tm2 = (1 - ICC.2 - ICC.3)*(1-R2.1) / ( (Tbar * (1-Tbar) * K * nbar) )
        J <- Q * (tm1 + tm2)
    } else if (d_m == 'd2.2_m2rc')
    {
        numr <- nbar * ICC.2 * (1 - R2.2) +
          (1 - ICC.2) * (1 - R2.1)
        denom <- Tbar * (1 - Tbar) * nbar
        J <- (MT/MDES)^2 * numr/denom
    } else if (d_m == 'd3.3_m3rc2rc')
    {
        numr <- nbar * ICC.2 * (1 - R2.2) +
          (1 - ICC.2 - ICC.3) * (1 - R2.1)
        denom <- nbar * ( Tbar * (1 - Tbar) * K *
          (MDES/MT)^2 - ICC.3 * (1 - R2.3) )
        J <- numr/denom
    } else if (d_m == 'd3.2_m3ff2rc' || d_m == 'd3.2_m3fc2rc' )
    {
        numr <- nbar * ICC.2 * (1 - R2.2) +
          (1 - ICC.2 - ICC.3) * (1 - R2.1)
        denom <- nbar * Tbar * (1 - Tbar) *
          K * (MDES/MT)^2
        J <- numr/denom
    } else if (d_m == 'd3.2_m3rr2rc')
    {
        numr <- nbar * ICC.2 * (1 - R2.2) +
            (1 - ICC.2 - ICC.3) * (1 - R2.1)
        denom <- nbar * Tbar * (1 - Tbar) *
            ( K * (MDES/MT)^2 - ICC.3 * omega.3 )
        J <- numr/denom
    } else
    {
        stop(paste('d_m not implemented:', d_m))
    }
    return( J )
}

#' Calculates K, the number of districts
#'
#' @param d_m a single RCT d_m (see list/naming convention)
#' @param MT multiplier
#' @param MDES scalar; the MDES value for all outcomes
#' @param J scalar; the number of schools
#' @param nbar scalar; the harmonic mean of the number of units per school
#' @param Tbar scalar; the proportion of samples that
#' are assigned to the treatment
#' @param R2.1 scalar, or vector of length M;
#' percent of variation explained by Level 1 covariates for each outcome
#' @param R2.2 scalar, or vector of length M;
#' percent of variation explained by Level 2 covariates for each outcome
#' @param R2.3 scalar, or vector of length M;
#' percent of variation explained by Level 3 covariates for each outcome
#' @param ICC.2 scalar; school intraclass correlation
#' @param ICC.3 scalar; district intraclass correlation
#' @param omega.2 scalar; ratio of school effect
#' size variability to random effects
#'   variability
#' @param omega.3 scalar; ratio of district
#' effect size variability to random effects
#'   variability
#'
#' @return K, the number of districts
#' @keywords internal
calc_K <- function(d_m, MT = 2.8, MDES, J, nbar, Tbar,
                   R2.1, R2.2, R2.3,
                   ICC.2, ICC.3,
                   omega.2, omega.3) {

    K <- NA
    if (d_m == 'd3.1_m3rr2rr')
    {
        K <- (MT/MDES)^2 *
            ( (ICC.3 * omega.3) +
                  (ICC.2 * omega.2) / J +
                  ((1 - ICC.2 - ICC.3) * (1 - R2.1))/(Tbar * (1 - Tbar) * J * nbar) )


        Q = (MT/MDES)^2
        tm1 = (ICC.2 * omega.2) / K
        tm2 = (1 - ICC.2 - ICC.3)*(1-R2.1) / ( (Tbar * (1-Tbar) * K * nbar) )
        J <- Q * (tm1 + tm2)
    } else if (d_m == 'd3.1_m3ff2rr')
    {
        Q = (MT/MDES)^2
        tm1 = (ICC.2 * omega.2) / J
        tm2 = (1 - ICC.2 - ICC.3)*(1-R2.1) / ( (Tbar * (1-Tbar) * J * nbar) )
        K <- Q * (tm1 + tm2)
    } else if (d_m == 'd3.3_m3rc2rc')
    {
        K <- (MT/MDES)^2 *
          ( (ICC.3 * (1 - R2.3)) / (Tbar * (1 - Tbar)) +
          (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
          ((1 - ICC.2 - ICC.3)*(1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
    } else if (d_m == 'd3.2_m3ff2rc' || d_m == 'd3.2_m3fc2rc' )
    {
        K <- (MT/MDES)^2 *
          ( (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
          ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
    } else if (d_m == 'd3.2_m3rr2rc')
    {
        K <- (MT/MDES)^2 *
          ( (ICC.3 * omega.3) +
          (ICC.2 * (1 - R2.2)) / (Tbar * (1 - Tbar) * J) +
          ((1 - ICC.2 - ICC.3) * (1 - R2.1)) / (Tbar * (1 - Tbar) * J * nbar) )
    } else
    {
        stop(paste('d_m not implemented:', d_m))
    }
    return(K)
}

#### Parameter and call validation code ####

make_MDES_vector <- function(MDES, M,
                             numZero = NULL,
                             verbose = TRUE) {

    if ( !is.null(numZero) ) {
        if ( numZero >= M )
        {
            stop('numZero (or propZero*M) cannot be greater than or equal to M')
        }
        if ( ( length(MDES) > 1 ) && ( numZero + length(MDES) != M ) )
        {
            stop(
             "Please provide an MDES vector + numZero (or propZero)
             that identify exactly M outcomes.\n
             Example: MDES = c(0.1, 0.1), numZero = 3 or
             propZero = 3/5 with M = 5\n
             Assumed MDES vector = c(0.1, 0.1, 0, 0, 0)"
            )
        }
        if ( length(MDES) == 1 ) {
            MDES <- c(rep( MDES, M - numZero), rep(0, numZero) )
        } else {
            MDES <- c(MDES, rep(0, numZero))
        }
        if ( verbose ) {
            message('Assumed full MDES vector:', 'c(',
                    paste(MDES, collapse = ', '), ')')
        }
    }

    if (length(MDES) != M)
    {
        if ( length(MDES) == 1 ) {
            MDES <- rep( MDES, M )
        } else {
            stop(paste('Please provide a vector of MDES values of length 1 or',
                       'M. Current vector: ',
                       paste0( MDES, collapse = ", " ), 'M =', M))
        }
    }

    return( MDES )
}


validate_MTP <- function(
    MTP, power.call, mdes.call, ss.call, M, pdef, multi.MTP.ok = FALSE
) {

    if ( !multi.MTP.ok && length( MTP ) > 1 )
    {
        stop( 'Please provide only a single MTP procedure.' )
    }

    if ( !is.null( MTP ) && any( MTP == "raw" ) ) {
        MTP[ MTP == "raw" ] <- "None"
    }

    if (M == 1)
    {
        if ( !is.null( MTP ) && (MTP != "None" ) )
        {
            warning("Multiple testing corrections are not needed when M = 1.")
        }
        MTP <- "None"
    } else {
        if (is.null(MTP))
        {
            stop('Please provide a multiple test procedure (MTP).')
        } else if ( (mdes.call || ss.call) &&
                   any( MTP == 'None' ) && !pdef$indiv )
        {
            stop(paste('For all minimum or complete power specifications,',
                       'you must provide a MTP.'))
        } else if ( length( MTP ) == 1 && MTP == 'None' )
        {
            warning('Proceeding with multiple outcomes and no MTP.')
        }
    }

    chk <- MTP %in% pump_info()$Adjustment$Method

    if ( !all( chk ) ) {
        if ( length( MTP ) > 1 ) {
            msg <- sprintf( 'You have at least one invalid MTP: %s',
                             paste( "'", MTP[!chk], "'", sep = "",
                               collapse = ", " ) )
        } else {
            msg <- sprintf( '"%s" is an invalid MTP.', MTP )
        }
        stop( msg )
    }

    return( MTP )
}



#' Validate d_m string
#'
#' Ensure d_m is a supported pair of design and model.
#' If d_m is just a design, select a default model.
#' Convert PowerUp! names to our naming system as needed.
#'
#' @return Full d_m string that will be found in `pump_info()`
#' @keywords internal
validate_d_m <- function(d_m) {
    # allow either supported d_m names or PowerUp! equivalents
    info <- pump_info()
    if ( !(d_m %in% info$Context$d_m) ) {
        if (d_m %in% info$Context$PowerUp) {
            d_m <- info$Context$d_m[info$Context$PowerUp == d_m]
        } else {
            dm <- parse_design(d_m)
            if ( is.null( dm ) ) {
                stop( glue::glue( '{d_m} is an invalid d_m.') )
            } else {
                match_index <- which(info$Context$Design == dm$design)
                if ( length( match_index ) == 0 ) {
                    stop( glue::glue( '{d_m} is an invalid d_m.') )
                } else {
                    options <- paste0( info$Context$d_m[ match_index ],
                                       collapse = ", " )
                    d_m <- info$Context$d_m[[match_index[[1]]]]
                    warning(glue::glue(paste("Selecting design and model",
                                             "{d_m} as default for design",
                                             "from following options: {options}")),
                        call. = FALSE )
                }
            }
        }
    }

    return( d_m )
}


#' Validates user inputs
#'
#' This functions takes in a list of user inputs. Depending on the inputs,
#' it produces errors or warnings, and at times modifies inputs if necessary.
#'
#' @param d_m a single RCT d_m (see list/naming convention)
#' @param params.list a list of parameters input by a user
#' @param power.call flag for power estimation
#' @param ss.call flag for sample size estimation
#' @param mdes.call flag for MDES estimation
#' @param verbose whether to print out warnings
#' @param multi.MTP.ok whether validation allows for multiple MTPs
#'
#' @return params.list
#' @keywords internal
validate_inputs <- function(d_m, params.list,
                            power.call = FALSE,
                            mdes.call = FALSE,
                            ss.call = FALSE,
                            verbose = TRUE,
                            multi.MTP.ok = FALSE)
{

    verbose_message <- function(msg)
    {
        if (verbose)
        {
            warning( msg )
        }
    }

    #-------------------------------------------------------#
    # basic checks of inputs
    #-------------------------------------------------------#

    d_m = validate_d_m( d_m )

    par.d_m <- parse_d_m(d_m)
    pdef <- parse_power_definition(
        params.list$power.definition, params.list$M
    )

    params.list$MTP <- validate_MTP( MTP = params.list$MTP,
                                    power.call = power.call,
                                    mdes.call = mdes.call,
                                    ss.call = ss.call,
                                    M = params.list$M,
                                    pdef = pdef,
                                    multi.MTP.ok = multi.MTP.ok )

    #-------------------------------------------------------#
    # Westfall-Young
    #-------------------------------------------------------#

    if ( ( any(params.list$MTP == "WY-SD") ||
           any(params.list$MTP == "WY-SS") ) &&
         params.list$B < 1000 )
    {
        warning(paste("For the step-down Westfall-Young procedure,",
                      "it is recommended that sample (B) be at least",
                      "1000. Current B:", params.list$B))
    }

    #-------------------------------------------------------#
    # MDES
    #-------------------------------------------------------#

    if ( !is.null(params.list$propZero) ) {
        if (!is.null(params.list$numZero) ) {
            stop( "Only one of numZero and propZero can be non-null" )
        }
        if ( params.list$propZero < 0 || params.list$propZero >= 1 ) {
            stop( "propZero must be in [0,1) range" )
        }
        params.list$numZero = round( params.list$M * params.list$propZero )
        params.list$propZero = NULL
    }

    if ( mdes.call ) {
        if ( !is.null( params.list$MDES ) ) {
            stop( "You cannot provide MDES to pump_mdes()" )
        }

        if ( !is.null( params.list$numZero ) &&
             params.list$numZero >= params.list$M ) {
            stop( sprintf(paste("You cannot specify %s zeros via numZero",
                                "or propZero with only %s outcomes"),
                 params.list$numZero, params.list$M ) )
        }

    } else {

        params.list$MDES <- make_MDES_vector(
            params.list$MDES,
            params.list$M, params.list$numZero,
            verbose = verbose )
    }

    #---------------------------------------------------------------#
    # enforce scalar parameters
    #---------------------------------------------------------------#
    if (!is.null(params.list$numCovar.1) && length(params.list$numCovar.1) != 1)
    {
        stop("numCovar.1: Please provide a scalar.")
    }
    if (!is.null(params.list$numCovar.2) && length(params.list$numCovar.2) != 1)
    {
        stop("numCovar.2: Please provide a scalar.")
    }
    if (!is.null(params.list$numCovar.3) && length(params.list$numCovar.3) != 1)
    {
        stop("numCovar.3: Please provide a scalar.")
    }


    #---------------------------------------------------------------#
    # convert all params from scalar to vector, if they are non-null
    #---------------------------------------------------------------#

    if (!(length(params.list$R2.1) %in% c(1, params.list$M)))
    {
        stop("R2.1: Please provide a scalar parameter or
             a vector of length M.")
    }
    if (length(params.list$R2.1) == 1)
    {
        params.list$R2.1 <- rep(params.list$R2.1, params.list$M)
    }

    if (!(length(params.list$R2.2) %in% c(0, 1, params.list$M)))
    {
        stop("R2.2: Please provide a scalar parameter or a vector of length M.")
    }
    if (length(params.list$R2.2) == 1)
    {
        params.list$R2.2 <- rep(params.list$R2.2, params.list$M)
    }

    if (!(length(params.list$R2.3) %in% c(0, 1, params.list$M)))
    {
        stop("R2.3: Please provide a scalar parameter or a vector of length M.")
    }
    if (length(params.list$R2.3) == 1)
    {
        params.list$R2.3 <- rep(params.list$R2.3, params.list$M)
    }

    if (!(length(params.list$ICC.2) %in% c(0, 1, params.list$M)))
    {
        stop("ICC.2: Please provide a scalar parameter or a vector of length M.")
    }
    if (length(params.list$ICC.2) == 1)
    {
        params.list$ICC.2 <- rep(params.list$ICC.2, params.list$M)
    }

    if (!(length(params.list$ICC.3) %in% c(0, 1, params.list$M)))
    {
        stop("ICC.3: Please provide a scalar parameter or a vector of length M.")
    }
    if (length(params.list$ICC.3) == 1)
    {
        params.list$ICC.3 <- rep(params.list$ICC.3, params.list$M)
    }

    if (!(length(params.list$omega.2) %in% c(0, 1, params.list$M)))
    {
        stop("omega.2: Please provide a scalar parameter or a vector of length M.")
    }
    if (length(params.list$omega.2) == 1)
    {
        params.list$omega.2 <- rep(params.list$omega.2, params.list$M)
    }

    if (!(length(params.list$omega.3) %in% c(0, 1, params.list$M)))
    {
        stop("omega.3: Please provide a scalar parameter or a vector of length M.")
    }
    if (length(params.list$omega.3) == 1)
    {
        params.list$omega.3 <- rep(params.list$omega.3, params.list$M)
    }

    #-------------------------------------------------------#
    # Basic checks of data parameters
    #-------------------------------------------------------#

    if ( ( !is.null( params.list$K ) && params.list$K <= 0) |
        ( !is.null( params.list$J ) && params.list$J <= 0) |
        ( !is.null( params.list$nbar ) && params.list$nbar <= 0) ) {
        stop('Provided values of J, K, and/or nbar need to be positive.')
    }

    if ( params.list$numCovar.1 < 0 |
        ( !is.null( params.list$numCovar.2 )  && params.list$numCovar.2 < 0  ) |
        ( !is.null( params.list$numCovar.3 ) && params.list$numCovar.3 < 0 ) )
    {
        stop('Please provide non-negative values of your num.Covar parameters')
    }

    if (params.list$Tbar >= 1 | params.list$Tbar <= 0)
    {
        stop('Please provide Tbar as a probability strictly between 0 and 1')
    }

    if (params.list$alpha > 1 | params.list$alpha < 0)
    {
        stop('Please provide alpha as a probability between 0 and 1')
    }

    if (any(params.list$R2.1 > 1) | any(params.list$R2.1 < 0) |
       any(params.list$R2.2 > 1) | any(params.list$R2.2 < 0) |
       any(params.list$R2.3 > 1) | any(params.list$R2.3 < 0))
    {
        stop('Please provide R2 as a probability between 0 and 1')
    }

    if (any(params.list$omega.2 < 0) | (!is.null(params.list$omega.3) &&
       any(params.list$omega.3 < 0)))
    {
        stop('Please provide a non-negative value of Omega')
    }

    # ICC
    if (!is.null(params.list$ICC.2) && !is.null(params.list$ICC.3) &&
       any(params.list$ICC.2 + params.list$ICC.3 > 1))
    {
      stop('ICC.2 + ICC.3 must be <= 1')
    }

    # ICC
    if (!is.null(params.list$ICC.2) && !is.null(params.list$ICC.3) &&
       any(params.list$ICC.2 + params.list$ICC.3 == 1))
    {
        warning('ICC.2 + ICC.3 = 1, leaving no variation at level 1')
    }

    if ( is.null( params.list$rho ) ) {
        if ( params.list$M > 1 ) {
            stop('Please provide rho as a correlation between -1 and 1')
        } else {
            params.list$rho = 0
        }
    } else {
        # invalid rho
        if (params.list$rho > 1 | params.list$rho < -1)
        {
          stop('Please provide rho as a correlation between -1 and 1')
        }
    }

    #-------------------------------------------------------#
    # check for inconsistent user inputs
    #-------------------------------------------------------#

    # one level models
    if ( par.d_m$levels == 1 ) {
      if ((!is.null( params.list$J ) && params.list$J != 1 ) ||
          (!is.null( params.list$K ) && params.list$K != 1 ) ||
          (!is.null( params.list$numCovar.2 ) && params.list$numCovar.2 > 0 ) ||
          (!is.null( params.list$R2.2 ) && any(params.list$R2.2 > 0 ) ) ||
          (!is.null( params.list$ICC.2 ) && any(params.list$ICC.2 > 0 ) ) ||
          (!is.null( params.list$omega.2 ) && any(params.list$omega.2 > 0 ) ) )

        warning(paste('The following parameters are not valid for one-level',
                      'designs, and will be ignored:\n',
                      'J, K, numCovar.2, R2.2, ICC.2, omega.2'))
      params.list$J <- NULL
      params.list$R2.2 <- NULL
      params.list$ICC.2 <- NULL
      params.list$omega.2 <- NULL
    }

    # check for three-level parameters when working with
    # one or two level models
    if ( par.d_m$levels <= 2 )
    {
      if ( ( !is.null(params.list$K) && params.list$K > 1 ) |
          ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
          ( !is.null(params.list$R2.3)) && any( params.list$R2.3 > 0 ) |
          ( !is.null(params.list$ICC.3)) && any( params.list$ICC.3 > 0 ) |
          ( !is.null(params.list$omega.3) && any(params.list$omega.3 > 0 ) ) )
      {
        warning(paste('The following parameters are only valid for three-level',
                      'designs, and will be ignored:\n',
                      'K, numCovar.3, R2.3, ICC.3, omega.3'))
        params.list$K <- NULL
        params.list$R2.3 <- NULL
        params.list$ICC.3 <- NULL
        params.list$omega.3 <- NULL
      }
    }

    # two level models and above
    if ( par.d_m$levels >= 2 )
    {
      if ( params.list$J == 1 )
      {
            warning('Running a model with >= 2 levels with J = 1')
      }
      if ( any(params.list$ICC.2 == 0 ) )
      {
            verbose_message('Running a model with >= 2 levels with ICC.2 = 0')
      }
      if ( par.d_m$rand_level == 2 && any( params.list$R2.2 == 0 ) )
      {
            verbose_message('Assuming R2.2 = 0')
      }
      if ( par.d_m$FE.2 &
          (( !is.null(params.list$numCovar.2) && params.list$numCovar.2 > 0 ) |
           ( !is.null(params.list$R2.2) && any( params.list$R2.2 > 0 ) ) ))
      {
        warning(paste('The following parameters are not valid for fixed effect',
                      'designs, and will be ignored:\n',
                      'numCovar.2, R2.2'))
            params.list$R2.2 <- NULL
      }
      if ( par.d_m$model2.p[2] == 'r' && any( params.list$omega.2 == 0 ) )
      {
            verbose_message('Fitting random impact model when omega.2 = 0')
      }
      # constant treatment effects models: level 2
      if ( par.d_m$model2.p[2] == 'c' )
      {
        if (any(params.list$omega.2 > 0))
        {
            verbose_message(paste('omega.2 will be ignored for constant',
                                  'treatment effects models.\n',
                                  'Ignoring input omega.2 value'))
        }
      }
    }


    # three level models
    if ( par.d_m$levels == 3 )
    {
      if (is.null(params.list$K) || params.list$K < 1 )
      {
        stop(paste('You must specify K, with K >= 1 (number of units',
                   'at level 3) for three-level designs' ))
      }
      if ( params.list$K == 1 )
      {
         warning('Running a 3-level model with K = 1')
      }
      if ( any( params.list$ICC.3 == 0 ) )
      {
          verbose_message('Running a 3-level model with ICC.3 = 0')
      }
      if ( par.d_m$rand_level == 3 & any(params.list$R2.3 == 0 ) )
      {
          verbose_message('Assuming R2.3 = 0')
      }
      if ( par.d_m$model3.p[2] == 'r' & any(params.list$omega.3 == 0 ) )
      {
          verbose_message('Assuming omega.3 = 0')
      }
      if ( par.d_m$FE.3 &
         (( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
          ( !is.null(params.list$R2.3) && any( params.list$R2.3 > 0 ) ) ))
      {
        warning(paste('The following parameters are not valid for fixed',
                      'effect designs, and will be ignored:\n',
                      'numCovar.3, R2.3'))
        params.list$R2.3 <- NULL
      }

      # constant treatment effects models: level 3
      if ( par.d_m$model3.p[[2]] == 'c' )
      {
        if (any(params.list$omega.3 > 0))
        {
          warning(paste('Omega is assumed to be 0 for constant treatment',
                        'effects models.\n',
                        'Ignoring input omega.3 value'))
          params.list$omega.3 <- NULL
        }
      }

    }

    # number covariates
    if (!is.null( params.list$R2.1) && any(params.list$R2.1 != 0) &&
                 params.list$numCovar.1 == 0)
    {
        warning(paste('If nonzero R2 (R2.1, level 1), at least one covariate',
                      'is assumed. Setting numCovar.1 = 1'))
        params.list$numCovar.1 <- 1
    }
    if (!is.null( params.list$R2.2) && any(params.list$R2.2 != 0) &&
                 params.list$numCovar.2 == 0)
    {
        warning(paste('If nonzero R2 (R2.2, level 2), at least one covariate',
                      'is assumed. Setting numCovar.2 = 1'))
        params.list$numCovar.2 <- 1
    }
    if (!is.null( params.list$R2.3) && any(params.list$R2.3 != 0) &&
                 params.list$numCovar.3 == 0)
    {
        warning(paste('If nonzero R2 (R2.3, level 3), at least one covariate',
                      'is assumed. Setting numCovar.3 = 1'))
        params.list$numCovar.3 <- 1
    }

    #-------------------------------------------------------#
    #  rho
    #-------------------------------------------------------#
    if (is.null(params.list$rho.matrix) && is.null(params.list$rho))
    {
        stop(sprintf(paste('Please provide either a %d x %d rho.matrix',
                           'or default scalar rho.'),
                     params.list$M, params.list$M ) )
    }

    if (!is.null(params.list$rho.matrix))
    {
        if (nrow(params.list$rho.matrix) != params.list$M |
           ncol(params.list$rho.matrix) != params.list$M)
        {
            stop(paste('Correlation matrix of invalid dimensions.',
                    'Please provide valid correlation matrix.'))
        }
        if (any(params.list$rho.matrix < -1) | any(params.list$rho.matrix > 1) )
        {
          stop('Please provide a rho matrix with values between -1 and 1.')
        }
    }

    params.list$propZero = NULL

    params.list$d_m = d_m

    return(params.list)

}


#' # Not yet implemented
#' # #' @param grid Flag of whether call is a grid call or non-grid call.
#'
#'
#' #' Check user inputs
#' #'
#' #' This functions takes in a list of user inputs and
#' #' checks them for validity,
#' #' producing a mix of errors or warnings.
#' #'
#' #' @param d_m a single RCT design (see list/naming convention)
#' #' @param call String denoting intended pump_power, pump_mdes, or pump_sample
#' #'   call.
#' #' @param verbose Keep and return smaller messages or not
#' #' @param multi.MTP.ok TRUE/FALSE on whether multiple MTP
#' #' is allowed for check.
#' #' @param ... The arguments to be passed to given call.
#' #'
#' #' @return list of two things: "ok", a TRUE/FALSE of whether this is a
#' #' valid call. "messages", a list of all messages, warnings, and notes
#' #' that would be generated by the input parameter choices.  In the case
#' #' of an error, the messages will have the error as
#' #' the first string in the list.
#' #' @import purrr
#' #' @export
#' check_pump_call = function( d_m,
#'                             call = c( "power", "mdes", "sample" ),
#'                             #grid = FALSE,
#'                             verbose = TRUE,
#'                             multi.MTP.ok = TRUE,
#'                             ... ) {
#'
#'     params = list( ... )
#'     call = match.arg(call)
#'     power.call = call == "power"
#'     mdes.call = call == "mdes"
#'     ss.call = call == "sample"
#'     params$d_m = d_m
#'
#'
#'     quiet_validate_inputs = silently( validate_inputs )
#'     cres = quiet_validate_inputs(
#'       d_m = d_m, params.list = params,
#'       power.call = power.call, mdes.call = mdes.call, ss.call = ss.call,
#'       verbose = verbose,
#'       multi.MTP.ok = multi.MTP.ok )
#'
#'     res = list()
#'     res$ok = is.null( cres$error )
#'
#'     messages = c()
#'     if ( !res$ok ) {
#'         messages = cres$error$message
#'     }
#'     messages = c( messages, cres$warnings, cres$messages )
#'     if ( cres$output != "" ) {
#'         messages = c( messages, cres$output )
#'     }
#'     res$messages = messages
#'     return( res )
#' }


