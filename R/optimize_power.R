#' #' Midpoint function
#' #'
#' #' Calculating the midpoint between the lower and upper bound by calculating
#' #' half the distance between the two and adding the lower bound to it. The
#' #' function is a helper function in determining the MDES that falls within
#' #' acceptable power range.
#' #'
#' #' @param lower lower bound
#' #' @param upper upper bound
#' #' @importFrom stats dist
#' #' @return returns midpoint value
#' midpoint <- function(lower, upper) {
#'   return(lower + dist(c(lower, upper))[[1]]/2)
#' }

#' Optimizes power to help in search for MDES or SS
#'
#' @param design a single RCT design (see list/naming convention)
#' @param search.type options: MDES, J, K (nbar not currently supported)
#' @param MTP a single multiple adjustment procedure of interest. Supported
#'   options: Bonferroni, BH, Holm, WY-SS, WY-SD
#' @param target.power Target power to arrive at
#' @param power.definition must be a valid power type outputted by power
#'   function, i.e. D1indiv, min1, etc.
#' @inheritParams pump_power
#'
#' @param cl cluster object to use for parallel processing
#' @param max.steps how many steps allowed before terminating
#' @param max.tnum maximum number of samples for a single step (other than the final check step).
#' @param final.tnum number of samples for final draw
#'
#' @return power
optimize_power <- function(design, search.type, MTP, target.power, power.definition, tol,
                           start.tnum, start.low, start.high,
                           MDES = NULL, J = NULL, K = 1, nbar = NULL,
                           M = M, Tbar = Tbar, alpha,
                           numCovar.1 = 0, numCovar.2 = 0, numCovar.3 = 0,
                           R2.1 = 0, R2.2 = 0, R2.3 = 0, ICC.2 = 0, ICC.3 = 0,
                           omega.2 = 0, omega.3 = 0, rho,
                           B = NULL, cl = NULL,
                           max.steps = 20, max.tnum = 2000, final.tnum = 4*max.tnum,
                           give.warnings = FALSE)
{
  
  # Helper function to call pump_power for our search and give back the power results
  # from the given set of parameters.
  power_check <- function( test_point, test_tnum ) {
    
    # Set K to passed parameter or test_point if we are searching on K
    # This is not below since ifelse() cannot allow a NULL value and K could be NULL for
    # 2 level designs.
    myK <- K
    if ( search.type == 'K' ) {
      myK <- test_point
    }
    
    if(search.type == 'mdes'){ MDES <- rep(test_point, M) }
    
    pt.power.results <- pump_power(
      design, MTP = MTP,
      MDES = MDES,
      nbar = ifelse(search.type == 'nbar', test_point, nbar),
      J = ifelse(search.type == 'J', test_point, J),
      K = myK,
      tnum = test_tnum,
      # fixed params
      M = M, Tbar = Tbar, alpha = alpha,
      numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
      R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3, ICC.2 = ICC.2, ICC.3 = ICC.3,
      rho = rho, omega.2 = omega.2, omega.3 = omega.3,
      B = B, cl = cl,
      validate.inputs = FALSE
    )
    
    return(pt.power.results)
  }
  
  power_check_df <- function( test_point, test_tnum ) {
    current.power.results <- power_check( test_point, test_tnum )
    
    current.power <- current.power.results[
      current.power.results$MTP == MTP, power.definition
    ]
    
    iter.results <- data.frame(
      step = step, 
      MTP = MTP, target.power = target.power,
      pt = test_point, w = test_tnum, power = current.power
    ) 
    return( iter.results )
  }
  
  
  gen_test_pts <- function(start.low, start.high, tnum)
  {
    # generate a series of points to try (on quadratic scale, especially relevant
    # for sample size)
    test.pts <- data.frame(
      step = 0,
      pt = seq(sqrt(start.low), sqrt(start.high), length.out = 5)^2,
      power = NA,
      w = tnum,
      MTP = MTP,
      target.power = target.power
    )
    
    # generate power for all the initial test points
    for(i in 1:nrow(test.pts))
    {
      pt.power.results <- power_check( test.pts$pt[i], start.tnum )
      if(is.na(pt.power.results$D1indiv[1]))
      {
        test.pts$power[i] <- 0
      } else
      {
        # Sanity check that we are getting power for a given metric.
        if(!(power.definition %in% colnames(pt.power.results)))
        {
          stop(paste0('Please provide a valid power definition. Provided definition: ', power.definition,
                      '. Available options: ', paste(colnames(pt.power.results), collapse = ', ')))
        }
        test.pts$power[i] <- pt.power.results[
          pt.power.results$MTP == MTP, power.definition
        ]
      }
    }
    
    return(test.pts)
  }
  
  # Ensure we have single MDES that is appropriate
  if ( search.type != "mdes" ) {
    stopifnot( !is.null( MDES ) )
    stopifnot( length(MDES) == 1 && MDES > 0 )
    MDES <- rep( MDES, M )
  }
  
  # Step 1: fit initial series of points to start search
  test.pts <- gen_test_pts(start.low, start.high, tnum = start.tnum)
  
  
  # Did we get NAs?  If so, currently crash (but should we impute 0 to keep
  # search going?)
  # TODO: test if we should keep going?
  # stopifnot( all( !is.na( test.pts$power ) ) )
  
  optimizer.warnings <- NULL
  quiet_find_best = purrr::quietly( find_best )
  
  # Based on initial grid, pick best guess for search.
  ct <- quiet_find_best( test.pts = test.pts, target.power = target.power, gamma = 1.5 )
  if ( !is.null( ct$warnings ) ) {
      optimizer.warnings <- c(optimizer.warnings, ct$warnings)
  }
  current.try = ct$result
  
  current.power <- 0
  current.tnum <- start.tnum
  step <- 0
  
  # Flag for if our next point to test is a valid (df > 0) point.  (Assume true
  # until we find otherwise.)
  current.try.ok <- TRUE
  
  bad_df_threshold = 0
  
  # Iteratively search by checking best point and then updating our curve.
  done = FALSE
  while( !done && step < max.steps )
  {
    step <- step + 1
    current.tnum <- pmin(max.tnum, round(current.tnum * 1.1))
    
    # what is smallest tested point?
    min_limit <- min( test.pts$pt )
    
    # the following code block is to handle possibly running out of degrees of
    # freedom if we are hugely overpowered.
    if ( current.try <= min_limit ) {
      # we are going lower than we have ever gone before.  Need to check if
      # degrees of freedom is still defined.
      
      current.try.ok = FALSE
      current.try = max( bad_df_threshold + 1, current.try )
      while( !current.try.ok && (min_limit - current.try) > 0.5 ) {
        
        myK = K
        if ( search.type == 'K' ) {
          myK = current.try
        }
        t.df <- calc.df(
          design = design,
          nbar = ifelse(search.type == 'nbar', current.try, nbar),
          J = ifelse(search.type == 'J', current.try, J),
          K = myK,
          numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3, validate=FALSE )
        
        if ( t.df > 1 ) {
          current.try.ok = TRUE
        } else {
          bad_df_threshold = max( bad_df_threshold, current.try )
          current.try = (current.try + min_limit) / 2
        }
      }
    } else {
      current.try.ok <- TRUE
    }
    
    if ( !current.try.ok ) {
      # We can't go lower than min_limit.  Need to see if it is indeed overpowered.
      current.try <- min_limit
    }
    
    iter.results <- power_check_df( current.try, current.tnum )
    
    # If we are close, check with more iterations and update our current step.
    if(abs(iter.results$power - target.power) < tol && current.tnum < max.tnum )
    {
      check.power.tnum <- pmin(10 * current.tnum, max.tnum - current.tnum)
      
      check.results <- power_check_df( current.try, check.power.tnum )
      
      avg_pow = weighted.mean( x = c(iter.results$power, check.results$power),
                               w = c(current.tnum, check.power.tnum) )
      
      # Overwrite results with our bonus step.
      iter.results$w <- current.tnum + check.power.tnum
      iter.results$power = avg_pow
    } # end if within tolerance
    
    # Record our step.
    current.power = iter.results$power
    test.pts <- dplyr::bind_rows(test.pts, iter.results)
    
    # If still good, or are stuck at minimum, go to a second full check to see
    # if we are winners!
    if( (abs(current.power - target.power) < tol) || 
        (!current.try.ok && current.power > target.power + tol) )
    {
      final.power.results <- power_check_df( current.try, final.tnum )
      
      test.pts <- dplyr::bind_rows(test.pts, final.power.results)
      current.power = final.power.results$power
    }
    
    if( (abs(current.power - target.power) < tol) || 
        (!current.try.ok && current.power > target.power + tol) ) {
      done = TRUE
    } else {
      
      ct <- quiet_find_best( test.pts = test.pts, target.power = target.power, gamma = 1.5 )
      if ( !is.null( ct$warnings ) ) {
        optimizer.warnings <- c(optimizer.warnings, ct$warnings)
      }
      current.try = ct$result
      
    }
  }
  
  # collect all warnings
  if(!is.null(optimizer.warnings) & give.warnings)
  {
    warning(unique(optimizer.warnings))
  }
  
  if ( !current.try.ok ) {
    warning( "Hit lower limit of what is allowed by degrees of freedom.  Likely overpowered." )
  }
  
  if( (step == max.steps) & abs(current.power - target.power) > tol) {
    msg <- "Reached maximum iterations without converging on estimate within tolerance.\n"
    msg <- paste(msg, "See sample size vignette for suggestions.")
    warning(msg)
    iter.results <- data.frame(
      step = step, MTP = MTP, target.power = target.power,
      pt = NA, power = NA, w = NA
    )
    test.pts <- dplyr::bind_rows(test.pts, iter.results )
    final.pts <- NULL
  } else
  {
    final.pts <- gen_test_pts(
      start.low = start.low,
      start.high = ceiling(test.pts$pt[nrow(test.pts)]),
      tnum = max.tnum
    )
    final.pts <- final.pts[, c('MTP', 'target.power', 'pt', 'w', 'power')]
  }
  test.pts = dplyr::relocate( test.pts, step, MTP, target.power, pt, w, power )
  return(list(test.pts = test.pts, final.pts = final.pts))
}






bounded_logistic_curve = function( x, params ) {
  if ( !is.null( names( params ) ) ) {
    beta0 <- params[["beta0"]]
    beta1 <- params[["beta1"]]
    pmin <- params[["pmin"]]
    pmax <- params[["pmax"]]
  } else {
    beta0 = params[[1]]
    beta1 = params[[2]]
    pmin = params[[3]]
    pmax = params[[4]]
  }
  pmin + (pmax - pmin)*plogis( beta0 + beta1*x )
}




d_bounded_logistic_curve = function( x, params ) {
  if ( !is.null( names( params ) ) ) {
    beta0 <- params[["beta0"]]
    beta1 <- params[["beta1"]]
    pmin <- params[["pmin"]]
    pmax <- params[["pmax"]]
  } else {
    beta0 = params[[1]]
    beta1 = params[[2]]
    pmin = params[[3]]
    pmax = params[[4]]
  }
  
  delta = pmax - pmin
  lin = -1 * (beta0 + beta1*x)
  e_lin = exp( lin )
  deriv = delta * ( 1/(1+e_lin)^2 ) * e_lin * beta1
  return( deriv )
}



find_crossover = function( target_power, params ) {
  beta0 <- params[["beta0"]]
  beta1 <- params[["beta1"]]
  pmin <- params[["pmin"]]
  pmax <- params[["pmax"]]
  
  delta = pmax - pmin
  xover = - ( log( delta / (target_power - pmin) - 1 ) + beta0 ) / beta1
  xover
}




#' Fit a bounded logistic curve
#' 
#' Curve is of form f(y) = pmin + (pmax-pmin) * logistic( beta0 + beta1*x )
#' 
#' (logistic as defined by plogis)
#' 
#' @return Vector of four estimated parameters for the logistic curve: beta0, beta1, pmin, pmax
fit_bounded_logistic = function( x, y, wt ) {
  
  # the log likelihood
  loglik <- function(par, y, x, wt) {
    p = bounded_logistic_curve( x, par )
    -sum( wt * (p - y)^2 )
  }
  
  # fit the model
  rs <- optim( c(0, 0.5 ,.1, .9), 
               loglik, 
               control=list(fnscale=-1),
               y=y, x=x, wt=wt, 
               lower=c(-Inf,0,0,0),
               upper=c(Inf,Inf,1,1),
               method="L-BFGS-B")
  names(rs$par) = c( "beta0", "beta1", "pmin", "pmax" )
  rs$par
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

find_best <- function(test.pts, target.power, gamma = 1.5)
{
  start.low <- min( test.pts$pt )
  start.high <- max( test.pts$pt )
  
  fit = fit_bounded_logistic( test.pts$pt, test.pts$power, sqrt( test.pts$w ) )
  
  if ( fit[["pmin"]] > target.power ) {
    # Min power is too high.  Reach lower.
    warning( "Minimum estimated power higher than target power" )
    return( start.low / gamma )
    
  } else if ( fit[["pmax"]] < target.power ) {
    # Max power is too low.  This could be a limitation or estimation error.
    warning( "Maximum estimated power lower than target power" )
    
    return( start.high * gamma )
  }
  
  #print( plot_power_search(test.pts) )
  
  # extract point where it crosses target power.
  cc = find_crossover( target.power, fit )
  return( cc )
}






find_best_semi_old <- function(test.pts, target.power, gamma = 1.5)
{
  start.low <- min( test.pts$pt )
  start.high <- max( test.pts$pt )
  
  # fit quadratic curve
  test.pts$sq_pt = sqrt( test.pts$pt )
  fit <- lm(
    power ~ 1 + sq_pt + I(sq_pt^2),
    data = test.pts,
    weights = w
  )
  
  # extract point where it crosses target power.
  # Our curve is now a x^2 + b x + (c - logit(target.power))
  # Using x = [ -b \pm sqrt( b^2 - 4a(c- logit(y))  ] / [2a]
  # first check if root exists
  cc <- rev( coef( fit ) )
  
  
  rt.check <- cc[2]^2 - 4 * cc[1] * (cc[3] - target.power)
  
  happy <- FALSE
  
  # We have a place where our quad line crosses target.power
  if ( rt.check > 0 ) {
    
    # calculate the two points to try (our two roots)
    try.pt <- ( -cc[2] + c(-1,1) * sqrt(rt.check) ) / 2 * cc[1]
    hits <- (start.low <= try.pt) & (try.pt <= start.high)
    
    if ( sum( hits ) == 1 ) {
      try.pt <- try.pt[hits]
      happy <- TRUE
    } else if ( sum( hits ) == 2 ) {
      # both roots in our range.  Probably flat.  Pick center
      try.pt = (try.pt[[1]] + try.pt[[2]]) / 2
      happy <- TRUE
    } else {
      happy <- FALSE
    }
  } else {
    warning( "No root in logistic model fit" )
  }
  
  if ( !happy ) {
    
    lin.mod <- lm( power ~ 1 + sq_pt, data = test.pts)
    cc <- rev( coef( lin.mod ) )
    try.pt <- ( (target.power - cc[[2]]) / cc[[1]] )^2
    
    if ( try.pt < start.low / gamma ) {
      try.pt <- start.low / gamma
    } else if ( try.pt > start.high * gamma ) {
      try.pt <- start.high * gamma
    }
  }
  
  return(unname(try.pt))
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

find_best_old <- function(test.pts, gamma = 1.5, target.power )
{
  # Get current range of search so far.
  start.low <- min( test.pts$pt )
  start.high <- max( test.pts$pt )
  
  # fit quadratic curve
  test.pts$sq_pt <- sqrt( test.pts$pt )
  quad.mod <- lm(
    power ~ 1 + sq_pt + I(sq_pt^2),
    weights = w,
    data = test.pts
  )
  
  # extract point where it crosses target power.
  # Our curve is now a x^2 + b x + (c - target.power)
  # Using x = [ -b \pm sqrt( b^2 - 4a(c-y) ) ] / [2a]
  # first check if root exists
  cc <- rev( coef( quad.mod ) )
  rt.check <- cc[2]^2 - 4 * cc[1] * (cc[3] - target.power)
  
  happy = FALSE
  if ( rt.check > 0 ) {
    # We have a place where our quad line crosses target.power
    
    # calculate the two points to try (our two roots)
    try.pt <- ( -cc[2] + c(-1,1) * sqrt(rt.check) ) / (2 * cc[1] )^2
    hits <- (start.low <= try.pt) & (try.pt <= start.high)
    
    if ( sum( hits ) == 1 ) {
      try.pt <- try.pt[hits]
      happy <- TRUE
    } else if ( sum( hits ) == 2 ) {
      # both roots in our range.  Probably flat.  Pick center
      try.pt = (try.pt[[1]] + try.pt[[2]]) / 2
      happy <- TRUE
    } else {
      happy <- FALSE
      # Go to linear
      
      # No roots in the original range.  Try extrapolation, if we can stay
      # within gamma of the original range
      try.pt <- sort( try.pt )
      if ( (try.pt[[2]] > start.low / gamma ) && ( try.pt[[1]] < start.high * gamma ) ) {
        if ( try.pt[[1]] < start.low / gamma ) {
          try.pt <- try.pt[[2]]
        } else if ( try.pt[[2]] > start.high * gamma ) {
          try.pt <- try.pt[[1]]
        } else {
          try.pt <- ifelse( sample(2,1) == 1, try.pt[[1]], try.pt[[2]] )
        }
        happy <- TRUE
      }
    }
  } else {
    warning( "No root in quadratic model fit" )
  }
  
  # If no roots in the original range, and quad extrapolation failed, try linear
  # extrapolation, but stay within gamma of the original range.
  if ( !happy ) {
    lin.mod <- lm( power ~ 1 + sq_pt, data = test.pts)
    cc <- rev( coef( lin.mod ) )
    try.pt <- ( (target.power - cc[[2]]) / cc[[1]] )^2
    
    if ( try.pt < start.low / gamma ) {
      try.pt <- start.low / gamma
    } else if ( try.pt > start.high * gamma ) {
      try.pt <- start.high * gamma
    }
  }
  
  return(try.pt)
}



#' Examine search path of the power search.
#'
#' This will give two plots about how the search narrowed down into the final
#' estimate.  Can be useful to gauge where convergence went poorly.
#'
#' @param pwr Result from the pump_sample or pump_mdes
#' 
#' @export
plot_power_search <- function( pwr ) {
  if ( is.pumpresult(pwr) ) {
    test.pts <- search_path(pwr)
  } else if ( is.data.frame(pwr) ) {
    test.pts = pwr
  } else {
    test.pts <- pwr$test.pts
  }
  
  require(gridExtra)
  require( ggplot2 )
  tp = dplyr::filter( test.pts, !is.na( power ) )
  fit = fit_bounded_logistic( tp$pt, tp$power, tp$w )
  
  lims = extendrange( r=range( test.pts$power, test.pts$target.power[[1]], na.rm=TRUE ), 0.15 )
  plot1 <-  ggplot( test.pts ) +
    geom_hline( yintercept = test.pts$target.power[1], col = "purple" ) +
    geom_point( aes( pt, power, size = w ), alpha = 0.5 ) + 
    theme_minimal() +
    geom_text( aes( pt, power, label = step ), vjust = "bottom", nudge_y = 0.01, size=3 ) +
    stat_function( col="red", fun = function(x) { bounded_logistic_curve( x, par=fit ) } ) +
    guides(colour="none", size="none") +
    coord_cartesian(ylim=lims )
  
  
  plot2 <-  ggplot( test.pts, aes(step, power, size = w) ) +
    geom_hline( yintercept = test.pts$target.power[1], col = "purple" ) +
    geom_point( alpha = 0.5 ) +
    scale_x_continuous( breaks=0:max(test.pts$step) ) +
    theme_minimal()+
    coord_cartesian(ylim=lims ) +
    guides(colour="none", size="none") 
  
  plot3 <-  ggplot( test.pts, aes(step, pt, size = w) ) +
    geom_point( alpha = 0.5 ) +
    scale_x_continuous( breaks=0:max(test.pts$step) ) +
    theme_minimal() +
    guides(colour="none", size="none") 
  
  grid.arrange(plot1, plot2, plot3, ncol=3) #+
  #title( "Search path for optimize_power" )
  
}

