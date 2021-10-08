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

    # Hack code since ifelse() cannot allow a NULL value and K could be NULL for
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

  gen_test_pts <- function(start.low, start.high, tnum)
  {
    # generate a series of points to try (on quadradic scale, especially relevant
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

  # Step 1: fit initial quadratic curve to start search
  test.pts <- gen_test_pts(start.low, start.high, tnum = start.tnum)


  # Did we get NAs?  If so, currently crash (but should we impute 0 to keep
  # search going?)
  # TODO: test if we should keep going?
  # stopifnot( all( !is.na( test.pts$power ) ) )

  optimizer.warnings <- NULL

  # Based on initial grid, pick best guess for search.
  current.try <- median( test.pts$pt )
  tryCatch(
    current.try <- find_best(
      test.pts = test.pts, target.power = target.power, gamma = 1.5
    ),
    warning = function(w) {
      optimizer.warnings <<- c(optimizer.warnings, w$message)
    }
  )
  current.power <- 0
  current.tnum <- start.tnum
  step <- 0

  # Flag if our next point to test is a valid (df > 0) point.  (Assume true
  # until we find otherwise.)
  current.try.ok <- TRUE

  # Iteratively search by checking best point and then updating our curve.
  while( (step < max.steps) & (abs( current.power - target.power ) > tol) )
  {
    step <- step + 1
    current.tnum <- pmin(max.tnum, round(current.tnum * 1.1))

    # what is smallest tested point?
    min_limit <- min( test.pts$pt )

    # This is a catch for running out of degrees of freedom if we are hugely overpowered.
    if ( current.try <= min_limit ) {
      # we are going lower than we have ever gone before.  Need to check if
      # degrees of freedom is still defined.

      if ( !current.try.ok ) {
        # we still want to go low.  We have hit a wall.
        break
      }

      current.try.ok = FALSE
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
          current.try = (current.try + min_limit) / 2
        }
      }
    } else {
      current.try.ok <- TRUE
    }

    if ( !current.try.ok ) {
      current.try <- min_limit
    }

    current.power.results <- power_check( current.try, current.tnum )

    current.power <- current.power.results[
      current.power.results$MTP == MTP, power.definition
    ]

    iter.results <- data.frame(
      step = step, pt = current.try, power = current.power, w = current.tnum,
      MTP = MTP, target.power = target.power
    )

    # If we are close, check with more iterations
    if(abs(current.power - target.power) < tol)
    {
      check.power.tnum <- pmin(10 * current.tnum, max.tnum)

      check.power.results <- power_check( current.try, check.power.tnum )
      check.power <- check.power.results[
        check.power.results$MTP == MTP, power.definition
      ]

      current.power <- (current.tnum*current.power + check.power.tnum*check.power)/(current.tnum+check.power.tnum)

      # Overwrite results with our bonus step.
      iter.results <- data.frame(
        step = step, pt = current.try, power = current.power,
        w = current.tnum + check.power.tnum,
        MTP = MTP, target.power = target.power
      )
    } # end if within tolerance

    # Record our step.
    test.pts <- dplyr::bind_rows(test.pts, iter.results)

    # If still good, go to a second full check to see if we are winners!
    if(abs(current.power - target.power) < tol)
    {
      final.power.results <- power_check( current.try, final.tnum )
      current.power <- final.power.results[
        final.power.results$MTP == MTP, power.definition
      ]

      iter.results <- data.frame(
        step = step, pt = current.try, power = current.power,
        w = final.tnum,
        MTP = MTP, target.power = target.power
      )
      test.pts <- dplyr::bind_rows(test.pts, iter.results)
    }

    tryCatch(
      current.try <- find_best(
        test.pts = test.pts, target.power = target.power, gamma = 1.5
      ),
      warning = function(w) {
        optimizer.warnings <<- c(optimizer.warnings, w$message)
      }
    )
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
      step = step, pt = NA, power = NA, w = NA,
      MTP = MTP, target.power = target.power
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
    final.pts <- final.pts[,c('pt', 'power', 'MTP', 'target.power')]
  }

  return(list(test.pts = test.pts, final.pts = final.pts))
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



#' Examine search path of the power search
#'
#' @param pwr Result from the pump_sample or pump_mdes
plot_power_search <- function( pwr, step_plot = FALSE ) {
  test.pts <- search_path(pwr)
  if ( !step_plot ) {
    ggplot( test.pts, aes( pt, power ) ) +
      geom_point( aes( size = w ), alpha = 0.5 ) + theme_minimal() +
      geom_hline( yintercept = test.pts$target.power[1], col = "red" ) +
      geom_text( aes( label = test.pts$step ), nudge_y = 0.01 )
  } else {
    ggplot( test.pts, aes(step, power, size = w) ) +
      geom_point() +
      geom_hline( yintercept = test.pts$target.power[1], col = "red" )
  }
}

