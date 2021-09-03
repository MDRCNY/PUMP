scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}

#' Validates user inputs
#'
#' This functions takes in a list of user inputs. Depending on the inputs,
#' it produces errors or warnings, and at times modifies inputs if necessary.
#'
#' @param design a single RCT design (see list/naming convention)
#' @param params.list a list of parameters input by a user
#'
#' @return params.list
#'
validate_inputs <- function( design, params.list,
                             mdes.call = FALSE,
                             single.MDES = FALSE)
{

  #-------------------------------------------------------#
  # basic checks of inputs
  #-------------------------------------------------------#

  if(!(design %in% c('d1.1_m2cc',
                     'd2.1_m2fc', 'd2.1_m2ff', 'd2.1_m2fr',
                     'd3.1_m3rr2rr', 'd2.2_m2rc', 'd3.3_m3rc2rc',
                     'd3.2_m3ff2rc', 'd3.2_m3rr2rc')))
  {
    stop('Invalid design.')
  }

  if(params.list$M == 1)
  {
    if ( !is.null( params.list$MTP ) && (params.list$MTP != "None" ) )
    {
      warning("Multiple testing corrections are not needed when M = 1.")
    }
    params.list$MTP <- "None"
  } else if(is.null(params.list$MTP))
  {
    stop('Please provide a multiple test procedure (MTP).')
  }

  if(length( params.list$MTP ) > 1)
  {
    stop( 'Please provide only a single MTP procedure.' )
  }

  if(!(params.list$MTP %in% c('None', 'Bonferroni', 'BH', 'Holm', 'WY-SS', 'WY-SD')))
  {
    stop('Invalid MTP.')
  }

  #-------------------------------------------------------#
  # MDES
  #-------------------------------------------------------#

  if ( mdes.call ) {
    if ( !is.null( params.list$MDES ) ) {
      stop( "You cannot provide MDES to pump_mdes()" )
    }
  } else
  {
    if(!is.null(params.list$numZero))
    {
      if( ( length(params.list$MDES) > 1 ) &&
          ( params.list$numZero + length(params.list$MDES) != params.list$M ) )
      {
        stop('Please provide an MDES vector + numZero that add up to M.\n
             Example: MDES = c(0.1, 0.1), numZero = 3, M = 5.\n
             Assumed MDES vector = c(0.1, 0.1, 0, 0, 0)')
      }
      if ( length(params.list$MDES) == 1 ) {
        params.list$MDES <- c(rep( params.list$MDES, params.list$M - params.list$numZero),
                              rep(0, params.list$numZero) )
      } else {
        params.list$MDES <- c(params.list$MDES, rep(0, params.list$numZero))
      }
      message('Assumed full MDES vector:', 'c(', paste(params.list$MDES, collapse = ', '), ')')
    }
    if ( single.MDES ) {
      if ( length(params.list$MDES) != 1 ) {
        stop( "Please provide a single MDES value. This function does not support vector MDES inputs." )
      }
    } else if(length(params.list$MDES) != params.list$M)
    {
      if ( length(params.list$MDES) == 1 ) {
        params.list$MDES <- rep( params.list$MDES, params.list$M )
        #warning('Assuming same MDES for all outcomes.  Specify full vector to remove this message.')
        #message('Assuming same MDES for all outcomes.  Specify full vector to remove this message.')
      } else {
        stop(paste('Please provide a vector of MDES values of length 1 or M. Current vector:',
                   MDES, 'M =', M))
      }
    }
  }

  #-------------------------------------------------------#
  # Basic checks of data parameters
  #-------------------------------------------------------#

  if( (!is.null( params.list$K ) && params.list$K <= 0) | ( !is.null( params.list$J) && params.list$J <= 0) | params.list$nbar <= 0)
  {
    stop('Pasted values of J, K, and/or nbar need to be positive.')
  }

  if(params.list$numCovar.1 < 0 | ( !is.null( params.list$numCovar.2 )  && params.list$numCovar.2 < 0  ) |
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

  if(params.list$omega.2 < 0 | (!is.null(params.list$omega.3) && params.list$omega.3 < 0))
  {
    stop('Please provide a non-negative value of Omega')
  }

  if(params.list$rho > 1 | params.list$rho < -1)
  {
    stop('Please provide rho as a correlation between -1 and 1')
  }

  if(!(length(params.list$R2.1) %in% c(1, params.list$M))  |
     !(length(params.list$R2.2) %in% c(1, params.list$M)) |
     !(length(params.list$R2.3) %in% c(1, params.list$M)) )
  {
    stop('Please provide R2 as a scalar or vector of length M.')
  }

  if(!(length(params.list$ICC.2) %in% c(1, params.list$M))  |
     !(length(params.list$ICC.3) %in% c(1, params.list$M)) )
  {
    stop('Please provide ICC as a scalar or vector of length M.')
  }

  if(!(length(params.list$omega.2) %in% c(1, params.list$M))  |
     !(length(params.list$omega.3) %in% c(1, params.list$M)) )
  {
    stop('Please provide omega as a scalar or vector of length M.')
  }
  #-------------------------------------------------------#
  # check for inconsistent user inputs
  #-------------------------------------------------------#

  if( design == 'd1.1_m2cc' ) {
    if ( !is.null( params.list$J ) && params.list$J != 1 ) {
      stop( "Can't have multiple units at 2nd level for the d1.1_m2cc design" )
    }

  #  if ( !is.null( params.list$numCovar.2 ) || !is.null( params.list$numCovar.3 ) || !is.null( params.list$R2.2 ) || !is.null( params.list$ICC.2 ) ) {
   #   stop( "Can't have numCovar.2, numCovar.3, R2.2, or ICC.2 for d1.1_m2cc design (single level design" )
  #  }

  }


  # two level models
  if(startsWith(design, 'd2') | startsWith(design, 'd1'))
  {
    if ( startsWith(design, 'd2') & params.list$J == 1 )
    {
      warning('Two level design with single unit at level 2')
    }

    if( ( !is.null(params.list$K) && params.list$K > 1 ) |
        ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
        ( !is.null(params.list$R2.3)) && any( params.list$R2.3 > 0 ) |
        ( !is.null(params.list$omega.3) && params.list$omega.3 > 0 ) )
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

  if(design == 'd3.2_m3ff2rc')
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

  # three level models
  if(startsWith(design, 'd3'))
  {
    if(is.null(params.list$K) || params.list$K <= 1 )
    {
      stop('You must specify K, with K > 1 (number of units at level 3) for three-level designs' )
    }
  }

  # three level models, continued.
  if(design %in% c('d3.1_m3rr2rr', 'd3.3_m3rc2rc', 'd3.2_m3rr2rc'))
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

  # constant treatment effects models: level 2
  if(design %in% c('d2.1_m2fc', 'd2.2_m2rc', 'd3.3_mrc2rc', 'd3.2_m3ff2rc'))
  {
    if(params.list$omega.2 > 0)
    {
      warning('Omega is assumed to be 0 for constant treatment effects models. Ignoring input omega.2 value')
      params.list$omega.2 <- 0
    }
  }
  # constant treatment effects models: level 3
  if(design %in% c('d3.3_mrc2rc'))
  {
    if(!is.null(params.list$omega.3) && params.list$omega.3 > 0)
    {
      warning('Omega is assumed to be 0 for constant treatment effects models. Ignoring input omega.3 value')
      params.list$omega.3 <- 0
    }
  }

  # specific 3 level models
  if(design %in% c('d3.1_m3rr2rr', 'd3.2_m3ff2rc', 'd3.2_m3rr2rc'))
  {
    if(is.null(params.list$omega.3))
    {
      stop('Omega.3 is required for this design.')
    }
  }
  if(design %in% c('d3.1_m3rr2rr', 'd3.2_m3rr2rc'))
  {
    if(is.null(params.list$ICC.3))
    {
      stop('ICC.3 is required for this design.')
    }
  }

  # number covariates
  if(!is.null( params.list$R2.1) && params.list$R2.1 != 0 && params.list$numCovar.1 == 0)
  {
    warning('If nonzero R2 (R2.1, level 1), at least one covariate is assumed. Setting numCovar.1 = 1')
    params.list$numCovar.1 <- 1
  }
  if(!is.null( params.list$R2.2) && params.list$R2.2 != 0 && params.list$numCovar.2 == 0)
  {
    warning('If nonzero R2 (R2.2, level 2), at least one covariate is assumed. Setting numCovar.2 = 1')
    params.list$numCovar.2 <- 1
  }
  if(!is.null( params.list$R2.3) && params.list$R2.3 != 0 && params.list$numCovar.3 == 0)
  {
    warning('If nonzero R2 (R2.3, level 3), at least one covariate is assumed. Setting numCovar.3 = 1')
    params.list$numCovar.3 <- 1
  }

  #-------------------------------------------------------#
  #  rho
  #-------------------------------------------------------#
  if(is.null(params.list$rho.matrix) & is.null(params.list$rho))
  {
    stop('Please provide either a rho matrix or default rho.')
  }

  if(!is.null(params.list$rho.matrix))
  {
    if(nrow(params.list$rho.matrix) != M | ncol(params.list$rho.matrix) != M)
    {
      stop('Correlation matrix of invalid dimensions. Please provide valid correlation matrix.')
    }
  }

  return(params.list)

}


#' Parse the power definition
#'
#' @param power.definition i.e. D1indiv, min1, complete
#' @param M number of outcomes
#' @return information about power type
parse_power_definition <- function( power.definition, M ) {
  powertype = list( min = FALSE,
                    complete = FALSE,
                    indiv = FALSE )

  if ( stringr::str_detect( power.definition, "min" ) ) {
    powertype$min = TRUE
    powertype$min_k = readr::parse_number( power.definition )
    stopifnot( is.numeric( powertype$min_k ) )
  } else if ( stringr::str_detect( power.definition, "complete" ) ) {
    powertype$min = TRUE
    powertype$complete = TRUE
    powertype$min_k = M
  } else if ( stringr::str_detect( power.definition, "indiv" ) ) {
    powertype$indiv = TRUE
    powertype$indiv_k = readr::parse_number( power.definition )
    stopifnot( is.numeric( powertype$indiv_k ) )
  }

  return( powertype )
}
