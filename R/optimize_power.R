


#' Optimizes power to help in search for MDES or SS
#'
#' @inheritParams pump_power
#' @inheritParams pump_sample
#' @param search.type type of optimization search being conducted. Options: K, J, nbar, mdes
#' @param start.low lower bound for optimization procedure
#' @param start.high upper bound for optimization procedure
#' @param give.warnings whether to return optimizer warnings
#' @param cl cluster object to use for parallel processing
#' @param grid.only TRUE means generate a grid from start.low to start.high, but
#'   do not do iterative search. (Useful for mapping out the power curve rather
#'   than identifying a point of particular power).
#' @param grid.size Number of points to check in initial search grid.  Grid will
#'   be spaced as a quadratic sequence (e.g., 0, 1, 4, 9, 16 for a 0-16 span).
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
                           give.warnings = FALSE,
                           grid.only = FALSE,
                           grid.size = 5 )
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

  # Bundle power_check results in the data frame
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

  # Generate grid of test points (no searching)
  gen_test_pts <- function(start.low, start.high, tnum, round=FALSE) {
    # generate a series of points to try (on quadratic scale, especially relevant
    # for sample size)
    pt = seq(sqrt(start.low), sqrt(start.high), length.out = grid.size)^2
    if ( round ) {
      pt <- unique( round(pt) )
    }
    test.pts <- data.frame(
      step = 0, pt = pt,
      w = tnum, MTP = MTP,
      target.power = target.power,
      power = NA
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
          stop(paste0(
            'Please provide a valid power definition. Provided definition: ', power.definition,
            '. Available options: ', paste(colnames(pt.power.results), collapse = ', ')))
        }
        test.pts$power[i] <- pt.power.results[
          pt.power.results$MTP == MTP, power.definition
        ]
      }
    }

    return(test.pts)
  }

  if ( search.type != "mdes" ) {
    start.low <- pmax( start.low, 1 )
  }


  # Step 1: fit initial series of points to start search
  test.pts <- gen_test_pts(start.low, start.high, tnum = start.tnum,
                           round = FALSE ) #search.type != "mdes" )

  if ( grid.only ) {
    return( test.pts )
  }

  optimizer.warnings <- NULL
  quiet_find_best <- purrr::quietly( find_best )

  # Based on initial grid, pick best guess for search.
  ct <- quiet_find_best( test.pts = test.pts, target.power = target.power, gamma = 1.5 )
  if ( !is.null( ct$warnings ) ) {
    optimizer.warnings <- c(optimizer.warnings, ct$warnings)
  }
  current.try <- ct$result$x
  current.try.dx <- ct$result$dx

  current.power <- 0
  current.tnum <- start.tnum
  step <- 0

  # Flag for if our next point to test is a valid (df > 0) point.  (Assume true
  # until we find otherwise.)
  current.try.ok <- TRUE

  bad_df_threshold <- 0

  # Iteratively search by checking best point and then updating our curve.
  done <- FALSE
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

      current.try.ok <- FALSE
      current.try <- max( bad_df_threshold + 1, current.try )
      while( !current.try.ok && (min_limit - current.try) > 0.5 ) {

        myK <- K
        if ( search.type == 'K' ) {
          myK <- current.try
        }
        t.df <- calc_df(
          design = design,
          nbar = ifelse(search.type == 'nbar', current.try, nbar),
          J = ifelse(search.type == 'J', current.try, J),
          K = myK,
          numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
          validate = FALSE
        )

        if ( t.df > 1 ) {
          current.try.ok <- TRUE
        } else {
          bad_df_threshold <- max( bad_df_threshold, current.try )
          current.try <- (current.try + min_limit) / 2
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
    iter.results$dx <- current.try.dx

    # If we are close, check with more iterations and update our current step.
    if(abs(iter.results$power - target.power) < tol && current.tnum < max.tnum )
    {
      check.power.tnum <- pmin(10 * current.tnum, max.tnum - current.tnum)

      check.results <- power_check_df( current.try, check.power.tnum )

      avg_pow <- stats::weighted.mean( x = c(iter.results$power, check.results$power),
                                       w = c(current.tnum, check.power.tnum) )

      # Overwrite results with our bonus step.
      iter.results$w <- current.tnum + check.power.tnum
      iter.results$power <- avg_pow
    } # end if within tolerance

    # Record our step.
    current.power <- iter.results$power
    test.pts <- dplyr::bind_rows(test.pts, iter.results)

    # If still good, or are stuck at minimum, go to a second full check to see
    # if we are winners!
    if( (abs(current.power - target.power) < tol) ||
        (!current.try.ok && current.power > target.power + tol) )
    {
      final.power.results <- power_check_df( current.try, final.tnum )
      final.power.results$dx <- current.try.dx
      test.pts <- dplyr::bind_rows(test.pts, final.power.results)
      current.power <- final.power.results$power

    }

    if( (abs(current.power - target.power) < tol) ||
        (!current.try.ok && current.power > target.power + tol) ) {
      done <- TRUE
    } else {

      ct <- quiet_find_best( test.pts = test.pts, target.power = target.power, gamma = 1.5 )
      if ( !is.null( ct$warnings ) ) {
        optimizer.warnings <- c(optimizer.warnings, ct$warnings)
      }
      current.try <- ct$result$x
      current.try.dx <- ct$result$dx
    }
  }

  # if ( current.try.dx < 0.0001 ) {
  #   warning( "Derivative flat; very possible estimate is off the mark.  Explore other (smaller) sample sizes" )
  # }

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
  }
  test.pts <- dplyr::relocate( test.pts, .data$step, .data$MTP, .data$target.power, .data$pt, .data$dx, .data$w, .data$power )
  return( test.pts = test.pts )
}


#' Calculate a power curve for sample size or mdes
#'
#' @description
#' For a grid of points based on a passed sample or mdes pumpresult, estimate
#' power.
#'
#' @param p pumpresult object
#' @param low Low end of grid
#' @param high High end of grid
#' @param grid.size Number of points in grid
#' @param tnum the number of test statistics (samples)
#'
#' @return List of powers for grid of points.
estimate_power_curve <- function( p, low = NULL, high = NULL,
                                  grid.size = 5, tnum = 2000 ) {
  stopifnot( is.pumpresult( p ) )
  stopifnot( pump_type(p) != "power" )

  pp <- params(p)
  pp$numZero <- NULL
  pp$rho.matrix <- NULL
  sp <- attr(p, "search.range" )
  test.pts <- search_path(p)
  if ( is.null( low ) ) {
    low <- sp[[1]]
  }
  if ( is.null( high ) ) {
    high <- sp[[2]] * 1.2
  }

  # for Bonferroni
  if(length(low) == 0)
  {
    low <- 0
  }
  if(length(high) == 0)
  {
    high <- p$`Sample.size`
  }

  search_type <- ifelse( pump_type(p) == "mdes",
                         "mdes",
                         attr(p, "sample.level" ) )
  final.pts <- do.call( optimize_power,
                        c( design = design(p),
                           pp, search.type = search_type,
                           start.low = low,
                           start.high = high,
                           start.tnum = tnum,
                           grid.only = TRUE,
                           grid.size = grid.size ) )

  return(final.pts)
}




bounded_logistic_curve <- function( x, params ) {
  if ( !is.null( names( params ) ) ) {
    beta0 <- params[["beta0"]]
    beta1 <- params[["beta1"]]
    pmin <- params[["pmin"]]
    pmax <- params[["pmax"]]
  } else {
    beta0 <- params[[1]]
    beta1 <- params[[2]]
    pmin <- params[[3]]
    pmax <- params[[4]]
  }
  pmin + (pmax - pmin)*stats::plogis( beta0 + beta1*x )
}




d_bounded_logistic_curve <- function( x, params ) {
  if ( !is.null( names( params ) ) ) {
    beta0 <- params[["beta0"]]
    beta1 <- params[["beta1"]]
    pmin <- params[["pmin"]]
    pmax <- params[["pmax"]]
  } else {
    beta0 <- params[[1]]
    beta1 <- params[[2]]
    pmin <- params[[3]]
    pmax <- params[[4]]
  }

  delta = pmax - pmin
  lin <- -1 * (beta0 + beta1*x)
  e_lin <- exp( lin )
  deriv <- delta * ( 1/(1+e_lin)^2 ) * e_lin * beta1
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
#' @param x The vector of covariate values of the logistics
#' @param y The proportion of 1s for a given value of x.  Same length as x.
#' @param wt The weight to place on a given x-y pair.  Same length as x, or
#'   scalar.
#'
#' @return Vector of four estimated parameters for the logistic curve: beta0,
#'   beta1, pmin, pmax
fit_bounded_logistic = function( x, y, wt ) {

  stopifnot( length(x) == length(y) )
  stopifnot( length(wt) == length(x) || length(wt) == 1 )

  # the log likelihood
  loglik <- function(par, y, x, wt) {
    p = bounded_logistic_curve( x, par )
    -sum( wt * (p - y)^2 )
  }

  # fit the model
  rs <- stats::optim( c(0, 0.5 ,.1, .9),
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
#' @param target.power goal power
#' @param gamma Number > 1. The amount we can extend our search range up (gamma)
#'   or down (1/gamma) if we want to search outside our current boundaries.  1
#'   means do not extend (which will potentially cause search to stall, not
#'   recommended).
#'
#' @return List of estimate of when curve reaches target.power, derivative of
#'   curve at that point, and parameters of the fit curve.

find_best <- function(test.pts, target.power, gamma = 1.5)
{
  stopifnot( gamma >= 1 )

  start.low <- min( test.pts$pt )
  start.high <- max( test.pts$pt )

  fit = fit_bounded_logistic( test.pts$pt, test.pts$power, sqrt( test.pts$w ) )

  if ( fit[["pmin"]] > target.power ) {
    # Min power is too high.  Reach lower.
    warning( "Minimum estimated power higher than target power" )
    cc =  start.low / gamma
  } else if ( fit[["pmax"]] < target.power ) {
    # Max power is too low.  This could be a limitation or estimation error.
    warning( "Maximum estimated power lower than target power" )

    cc = start.high * gamma
  } else {
    # extract point where it crosses target power.
    cc = find_crossover( target.power, fit )
  }

  return( list( x = cc, dx = d_bounded_logistic_curve(cc,fit), params=fit ) )
}






# find_best_semi_old <- function(test.pts, target.power, gamma = 1.5)
# {
#   start.low <- min( test.pts$pt )
#   start.high <- max( test.pts$pt )
#
#   # fit quadratic curve
#   test.pts$sq_pt = sqrt( test.pts$pt )
#   fit <- lm(
#     power ~ 1 + sq_pt + I(sq_pt^2),
#     data = test.pts,
#     weights = w
#   )
#
#   # extract point where it crosses target power.
#   # Our curve is now a x^2 + b x + (c - logit(target.power))
#   # Using x = [ -b \pm sqrt( b^2 - 4a(c- logit(y))  ] / [2a]
#   # first check if root exists
#   cc <- rev( coef( fit ) )
#
#
#   rt.check <- cc[2]^2 - 4 * cc[1] * (cc[3] - target.power)
#
#   happy <- FALSE
#
#   # We have a place where our quad line crosses target.power
#   if ( rt.check > 0 ) {
#
#     # calculate the two points to try (our two roots)
#     try.pt <- ( -cc[2] + c(-1,1) * sqrt(rt.check) ) / 2 * cc[1]
#     hits <- (start.low <= try.pt) & (try.pt <= start.high)
#
#     if ( sum( hits ) == 1 ) {
#       try.pt <- try.pt[hits]
#       happy <- TRUE
#     } else if ( sum( hits ) == 2 ) {
#       # both roots in our range.  Probably flat.  Pick center
#       try.pt = (try.pt[[1]] + try.pt[[2]]) / 2
#       happy <- TRUE
#     } else {
#       happy <- FALSE
#     }
#   } else {
#     warning( "No root in logistic model fit" )
#   }
#
#   if ( !happy ) {
#
#     lin.mod <- lm( power ~ 1 + sq_pt, data = test.pts)
#     cc <- rev( coef( lin.mod ) )
#     try.pt <- ( (target.power - cc[[2]]) / cc[[1]] )^2
#
#     if ( try.pt < start.low / gamma ) {
#       try.pt <- start.low / gamma
#     } else if ( try.pt > start.high * gamma ) {
#       try.pt <- start.high * gamma
#     }
#   }
#
#   return(unname(try.pt))
# }




# #' Extract roots from quadratic curve based on given evaluated points
# #'
# #' @param test.pts power evaluated at different points
# #' @param start.low lower bound
# #' @param start.high upper bound
# #' @param target.power goal power
# #' @param alternate alternate point to return if quadratic fit fails
# #'
# #' @return root of quadratic curve
#
# find_best_old <- function(test.pts, gamma = 1.5, target.power )
# {
#   # Get current range of search so far.
#   start.low <- min( test.pts$pt )
#   start.high <- max( test.pts$pt )
#
#   # fit quadratic curve
#   test.pts$sq_pt <- sqrt( test.pts$pt )
#   quad.mod <- lm(
#     power ~ 1 + sq_pt + I(sq_pt^2),
#     weights = w,
#     data = test.pts
#   )
#
#   # extract point where it crosses target power.
#   # Our curve is now a x^2 + b x + (c - target.power)
#   # Using x = [ -b \pm sqrt( b^2 - 4a(c-y) ) ] / [2a]
#   # first check if root exists
#   cc <- rev( coef( quad.mod ) )
#   rt.check <- cc[2]^2 - 4 * cc[1] * (cc[3] - target.power)
#
#   happy = FALSE
#   if ( rt.check > 0 ) {
#     # We have a place where our quad line crosses target.power
#
#     # calculate the two points to try (our two roots)
#     try.pt <- ( -cc[2] + c(-1,1) * sqrt(rt.check) ) / (2 * cc[1] )^2
#     hits <- (start.low <= try.pt) & (try.pt <= start.high)
#
#     if ( sum( hits ) == 1 ) {
#       try.pt <- try.pt[hits]
#       happy <- TRUE
#     } else if ( sum( hits ) == 2 ) {
#       # both roots in our range.  Probably flat.  Pick center
#       try.pt = (try.pt[[1]] + try.pt[[2]]) / 2
#       happy <- TRUE
#     } else {
#       happy <- FALSE
#       # Go to linear
#
#       # No roots in the original range.  Try extrapolation, if we can stay
#       # within gamma of the original range
#       try.pt <- sort( try.pt )
#       if ( (try.pt[[2]] > start.low / gamma ) && ( try.pt[[1]] < start.high * gamma ) ) {
#         if ( try.pt[[1]] < start.low / gamma ) {
#           try.pt <- try.pt[[2]]
#         } else if ( try.pt[[2]] > start.high * gamma ) {
#           try.pt <- try.pt[[1]]
#         } else {
#           try.pt <- ifelse( sample(2,1) == 1, try.pt[[1]], try.pt[[2]] )
#         }
#         happy <- TRUE
#       }
#     }
#   } else {
#     warning( "No root in quadratic model fit" )
#   }
#
#   # If no roots in the original range, and quad extrapolation failed, try linear
#   # extrapolation, but stay within gamma of the original range.
#   if ( !happy ) {
#     lin.mod <- lm( power ~ 1 + sq_pt, data = test.pts)
#     cc <- rev( coef( lin.mod ) )
#     try.pt <- ( (target.power - cc[[2]]) / cc[[1]] )^2
#
#     if ( try.pt < start.low / gamma ) {
#       try.pt <- start.low / gamma
#     } else if ( try.pt > start.high * gamma ) {
#       try.pt <- start.high * gamma
#     }
#   }
#
#   return(try.pt)
# }




#' Examine search path of the power search.
#'
#' This will give two plots about how the search narrowed down into the final
#' estimate.  Can be useful to gauge where convergence went poorly.
#'
#' @param pwr Result from the pump_sample or pump_mdes
#' @param plot.points flag; whether to plot individual points on curve
#'
#' @export
plot_power_curve <- function( pwr, plot.points = TRUE ) {
  if ( is.pumpresult( pwr ) ) {
    test.pts <- power_curve(pwr, all = TRUE )
    x_label <- pump_type(pwr)
  } else {
    stopifnot( is.data.frame(pwr) )
    test.pts <- pwr
    x_label <- "parameter"
  }

  tp <- dplyr::filter( test.pts, !is.na( .data$power ) )
  fit <- fit_bounded_logistic( tp$pt, tp$power, tp$w )

  xrng = range( test.pts$pt )
  lims <- grDevices::extendrange( r = range( test.pts$power, test.pts$target.power[[1]], na.rm = TRUE ), 0.15 )
  limsX <- grDevices::extendrange( r = xrng, 0.15 )

  plot1 <-  ggplot2::ggplot( test.pts ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[[1]], col = "purple" ) +
    ggplot2::theme_minimal() +
    ggplot2::stat_function( col="red", fun = function(x) { bounded_logistic_curve( x, params = fit ) } ) +
    ggplot2::guides(colour="none", size="none") +
    ggplot2::coord_cartesian( ylim=lims, xlim = limsX ) +
    ggplot2::labs( x = x_label, y = "power" )

  delrange = diff( xrng )
  if ( delrange >= 2 && delrange <= 10 ) {
    xpt = seq( floor( xrng[[1]] ), ceiling( xrng[[2]] ) )
    plot1 = plot1 +   ggplot2::scale_x_log10( breaks = xpt )
  } else {
    plot1 = plot1 +  ggplot2::scale_x_log10()
  }


  if ( plot.points ) {
    plot1 <- plot1 + ggplot2::geom_point( ggplot2::aes( .data$pt, .data$power, size = .data$w ), alpha = 0.5 )
  }

  return( plot1 )

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
    test.pts <- pwr
  } else {
    test.pts <- pwr$test.pts
  }

  if(is.null(test.pts))
  {
    stop('Algorithm converged in one iteration. No search path.')
  }

  tp <- dplyr::filter( test.pts, !is.na( .data$power ) )
  fit <- fit_bounded_logistic( tp$pt, tp$power, tp$w )

  lims <- grDevices::extendrange( r = range( test.pts$power, test.pts$target.power[[1]], na.rm=TRUE ), 0.15 )
  plot1 <-  ggplot2::ggplot( test.pts ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[1], col = "purple" ) +
    ggplot2::geom_point( ggplot2::aes( .data$pt, .data$power, size = .data$w ), alpha = 0.5 ) +
    ggplot2::theme_minimal() +
    ggplot2::geom_text( ggplot2::aes( .data$pt, .data$power, label = .data$step ), vjust = "bottom", nudge_y = 0.01, size=3 ) +
    ggplot2::stat_function( col = "red", fun = function(x) { bounded_logistic_curve( x, params =fit ) } ) +
    ggplot2::guides(colour = "none", size="none") +
    ggplot2::coord_cartesian(ylim = lims ) +
    ggplot2::scale_x_log10()


  plot2 <-  ggplot2::ggplot( test.pts, ggplot2::aes(.data$step, .data$power, size = .data$w) ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[1], col = "purple" ) +
    ggplot2::geom_point( alpha = 0.5 ) +
    ggplot2::scale_x_continuous( breaks=0:max(test.pts$step) ) +
    ggplot2::theme_minimal()+
    ggplot2::coord_cartesian(ylim=lims ) +
    ggplot2::guides(colour="none", size="none")

  plot3 <-  ggplot2::ggplot( test.pts, ggplot2::aes(.data$step, .data$pt, size = .data$w) ) +
    ggplot2::geom_point( alpha = 0.5 ) +
    ggplot2::scale_x_continuous( breaks=0:max(test.pts$step) ) +
    ggplot2::theme_minimal() +
    ggplot2::guides(colour="none", size="none")+
    ggplot2::scale_y_log10()

  gridExtra::grid.arrange(plot1, plot2, plot3, ncol=3) #+
  #title( "Search path for optimize_power" )

}
