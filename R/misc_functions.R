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
  
  #-------------------------------------------------------#
  # basic checks of inputs
  #-------------------------------------------------------#
  
  if(!(design %in% c('unblocked_i1_2c',
                     'blocked_i1_2c', 'blocked_i1_2f', 'blocked_i1_2r',
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
  
  #-------------------------------------------------------#
  # MDES
  #-------------------------------------------------------#
  
  if ( mdes.call ) {
    if ( !is.null( params.list$MDES ) ) {
      stop( "You cannot provide MDES to pump_mdes()" )
    }
  } else 
  {
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
  
  #-------------------------------------------------------#
  # Basic checks of data parameters
  #-------------------------------------------------------#
  
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
  
  #-------------------------------------------------------#
  # check for inconsistent user inputs
  #-------------------------------------------------------#
  
  if(M == 1)
  {
    warning("Multiple testing corrections are not needed when M = 1.")
    MTP <- "Bonferroni"
  }
  
  if(J == 1 & design != 'unblocked_i1_2c')
  {
    message('Assuming unblocked design')
    design <- 'unblocked_i1_2c'
  }
  
  # two level models
  if(design %in% c('blocked_i1_2c', 'blocked_i1_2f', 'blocked_i1_2r', 'simple_c2_2r'))
  {
    if( (!is.null(params.list$K) && params.list$K > 1 )|
        ( !is.null(params.list$numCovar.3) && params.list$numCovar.3 > 0 ) |
        ( !is.null(params.list$R2.3)) && any( params.list$R2.3 > 0 ) |
        ( !is.null(params.list$omega.3) && params.list$omega.3 > 0 ) )
    {
      warning('The following parameters are not valid for two-level designs, and will be ignored:\n
              K, numCovar.3, R2.3, ICC.3, omega.3')
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
      warning('The following parameters are not valid for fixed effect designs, and will be ignored:\n
              numCovar.3, R2.3')
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

  #-------------------------------------------------------#
  #  rho
  #-------------------------------------------------------#
  if(!is.null(rho.matrix) & !is.null(rho))
  {
    warning('Provided both rho and full rho matrix, using only rho.matrix')
    params.list$rho <- NULL
  }
  
  if(!is.null(rho.matrix))
  {
    if(nrow(rho.matrix) != M | ncol(rho.matrix) != M)
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
