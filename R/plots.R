

get_sample_tick_marks <- function( pt, breaks = 5, include.points = TRUE ) {
  
  pt <- round( pt )

  mn <- min( pt, na.rm = TRUE  )
  mx <- max( pt, na.rm = TRUE )
  
  sq <- round( seq( mn, mx, length.out = breaks ) )
  
  pts <- unique( c( sq, mn, mx ) )
  
  if ( include.points ) {
    pts <- union( pts, pt )
  }
  
  return( sort( pts ) )
}


#' @title Examine a power curve
#'
#' @description This will give a plot of power vs. 
#' MDES or sample size. It can be useful to
#' see how quickly power changes as a function of 
#' these design parameters. Can be useful to diagnose
#' relatively flat power curves, where power changes
#' little as a function of MDES or sample size, and can
#' also be useful to gauge where convergence went poorly.
#'
#' @param pwr Result from the pump_sample or pump_mdes (or data frame from,
#'   e.g., power_curve()).
#' @param plot.points flag; whether to plot individually tested points on curve
#' @param fit A four parameter bounded logistic curve
#' (if NULL will fit one to passed points).
#'
#' @inheritParams power_curve
#'
#' @export
#' 
#' @examples
#' mdes <- pump_mdes(d_m = "d2.1_m2fc", MTP = 'HO',
#'   power.definition = 'D1indiv', target.power = 0.7,
#'   J = 60, nbar = 50, M = 3, Tbar = 0.5, alpha = 0.05,
#'   numCovar.1 = 1, R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2)
#' plot_power_curve(mdes)
plot_power_curve <- function( pwr, plot.points = TRUE,
                              all = FALSE,
                              low = NULL, high = NULL,
                              grid.size = 5, tnum = 2000,
                              fit = NULL ) {
  
  
  if ( is.pumpresult( pwr ) ) {
    test.pts <- power_curve(pwr, low = low, high = high, 
                            grid.size = grid.size,
                            tnum = params(pwr)$tnum, all = all )
    x_label <- pump_type(pwr)
  } else {
    stopifnot( is.data.frame(pwr) )
    test.pts <- pwr
    x_label <- "parameter"
  }
  
  tp <- dplyr::filter( test.pts, !is.na( .data$power ) )
  
  if ( is.null(fit) ) {
    fit <- fit_bounded_logistic( tp$pt, tp$power, tp$w )
  }
  
  xrng <- range( test.pts$pt, na.rm = TRUE )
  limsX <- grDevices::extendrange( r = xrng, 0.15 )
  
  # sample size or MDES?  If dataframe passed, we don't know what we are
  # plotting so we turn to the scale of X to determine the lower bound.
  if ( limsX[1] <= 0 ) {
    if ( xrng[[2]] > 10 ) {
      limsX[1] <- 1
    } else {
      limsX[1] <- 0.0001
    }
  }
  
  plot1 <-  ggplot2::ggplot( test.pts ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[[1]],
                         col = "purple" ) +
    ggplot2::theme_minimal() +
    ggplot2::stat_function( col = "red",
      fun = function(x) { bounded_logistic_curve( x, params = fit ) } ) +
    ggplot2::guides(colour = "none", size="none") +
    ggplot2::coord_cartesian( ylim = c(0,1), xlim = limsX ) +
    ggplot2::labs( x = x_label, y = "power" )
  
  xpt <- get_sample_tick_marks(test.pts$pt, breaks = 5,
                               include.points = plot.points )
  
  delrange <- diff( xrng )
  if ( delrange > 50 ) {
    plot1 <- plot1 +  ggplot2::scale_x_log10( breaks=xpt )
  } else if ( delrange >= 2 && delrange <= 15 ) {
    # Tick marks for each sample size.
    xpt <- seq( floor( xrng[[1]] ), ceiling( xrng[[2]] ) )
    plot1 <- plot1 +
        ggplot2::scale_x_continuous( breaks = xpt )
  } else {
    plot1 <- plot1 +
        ggplot2::scale_x_continuous( breaks = xpt )
  }
  
  
  if ( plot.points ) {
    plot1 <- plot1 +
        ggplot2::geom_point(
         ggplot2::aes( .data$pt, .data$power, size = .data$w ), 
                       alpha = 0.5 )
  }
  
  return( plot1 )
  
}





#' @title Examine search path of a power search
#'
#' @description This will give triple-plots about 
#' how the search narrowed down into the
#' final estimate.  Can be useful to gauge where 
#' convergence went poorly.
#'
#' @param pwr Result from a pump_sample or pump_mdes call.
#' @param fit A fitted curve to the search.
#' @param target.line If non-NULL, add a reference line for the true power (if
#'   known, e.g., from a pump_power call). 
#'
#' @return a ggplot plot (a gridExtra arrangement of 3 plots, technically).
#' @export
#'
#' @examples
#'J <- pump_sample(d_m = "d2.1_m2fc",
#'    MTP = 'HO', power.definition = 'D1indiv',
#'    typesample = 'J', target.power = 0.6,
#'    nbar = 50, M = 3, MDES = 0.125,
#'    Tbar = 0.5, alpha = 0.05,
#'    numCovar.1 = 1, R2.1 = 0.1, ICC.2 = 0.05, 
#'    rho = 0.2, tnum = 2000)
#' plot_power_search(J)
#'
plot_power_search <- function( pwr, fit = NULL, target.line = NULL) {
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
  
  if ( is.null( fit ) ) {
    fit <- fit_bounded_logistic( tp$pt, tp$power, tp$w )
  }
  
  
  plot1 <-  plot_power_curve(test.pts, plot.points = TRUE )
  
  if ( !is.null( target.line ) ) {
    plot1 <- plot1 +
        ggplot2::geom_vline( xintercept = target.line, col = "purple" )
  }
  
  
  lims <- grDevices::extendrange( r = range( test.pts$power, 
                                             test.pts$target.power[[1]], 
                                             na.rm = TRUE ), 0.15 )
  plot2 <-  ggplot2::ggplot( test.pts,
    ggplot2::aes(.data$step, .data$power, size = .data$w) ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[1],
                         col = "purple" ) +
    ggplot2::geom_point( alpha = 0.5 ) +
    ggplot2::scale_x_continuous( breaks = 0:max(test.pts$step) ) +
    ggplot2::theme_minimal()+
    ggplot2::coord_cartesian(ylim = lims ) +
    ggplot2::guides(colour = "none", size = "none")
  
  plot3 <-  ggplot2::ggplot( test.pts,
    ggplot2::aes(.data$step, .data$pt, size = .data$w) ) +
    ggplot2::geom_point( alpha = 0.5 ) +
    ggplot2::scale_x_continuous( breaks=0:max(test.pts$step) ) +
    ggplot2::theme_minimal() +
    ggplot2::guides(colour="none", size="none")+
    ggplot2::scale_y_log10()
  
  if ( !is.null( target.line ) ) {
    plot3 <- plot3 +
        ggplot2::geom_hline( yintercept = target.line, col = "purple" )
  }
  
  gridExtra::grid.arrange(plot1, plot2, plot3, ncol=3)
  
}

#' @title Plot a single scenario pump object
#' 
#' @description Works on an object returned by pump_power(),
#' and visualizes different definitions of power across
#' MTPs. This function does not apply to pump_mdes()
#' or pump_sample() objects, as these functions only
#' return a single value.
#'
#' @param x pumpresult object
#' @param ... Additional parameters
#'
#' @return a ggplot plot.
#' @export
#'
#' @examples 
#' pp1 <- pump_power(d_m = "d2.2_m2rc", MTP = 'HO',
#'  nbar = 50, J = 20, M = 8, numZero = 5,
#'  MDES = 0.30, Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
#'  numCovar.1 = 1, numCovar.2 = 1, R2.1 = 0.1, R2.2 = 0.7, 
#'  ICC.2 = 0.05, rho = 0.2, tnum = 5000)
#'  
#' plot(pp1)

plot.pumpresult <- function( x, ... )
{
  stopifnot( is.pumpresult( x ) )
  
  if(pump_type(x) == 'power')
  {
    if(attr( x, "long.table" ))
    {
      x <- transpose_power_table(x)
    }
    
    plot.data <-
      x %>%
      dplyr::select_all() %>%
      dplyr::select(-.data$indiv.mean) %>%
      tidyr::pivot_longer(!.data$MTP,
                          names_to = "powerType", values_to = "power")
    
    # Creating power type as a factor for ordering on x axis
    M <- params(x)$M
    dPowers <- paste0("D",1:M,"indiv")
    minPowers <- paste0("min",1:M-1)
    complete <- "complete"
    plot.data$powerType <- factor(
      plot.data$powerType,
      levels = c(dPowers, minPowers, complete),
      ordered = TRUE
    )
    
    # remove missing rows
    plot.data <- plot.data %>%
      as.data.frame( plot.data ) %>%
      dplyr::filter( !is.na( .data$power ) )
    
    # single scenario plot
    ss.plot <- ggplot2::ggplot(
      data = plot.data,
      ggplot2::aes_string(x = "powerType",
                          y = "power",
                          shape = "MTP",
                          color = "MTP")) +
      ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(0.25)) +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::ggtitle(paste0("Adjusted power across
                               different definitions of power")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                                        face = "bold",
                                                        vjust = 1,
                                                        hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 10, angle = 45),
                     axis.text.y = ggplot2::element_text(size = 10),
                     axis.title  = ggplot2::element_text(size = 10)
      ) +
      ggplot2::labs(color = "", shape = "")
  } else if( pump_type(x) %in% c('mdes', 'sample') )
  {
    stop('plot() only works on pump_power() objects or grid objects,
         not pump_mdes() or pump_sample() objects.')
  } else
  {
    stop('Invalid pumpresult type.')
  }
  
  return(ss.plot)
}



#' Plot a grid pump power object
#'
#' @inheritParams plot.pumpgridresult
#' 
#' @return a ggplot object
#' 
#' @importFrom stringr str_detect
#' @keywords internal
plot.pumpgridresult.power <- function(
    x, power.definition = NULL, var.vary = NULL, ... 
) {
  
  M <- params(x)$M
  MTPs <- unique(c("None", params(x)$MTP))
  
  if(!attr( x, "long.table" )) {
    x <- transpose_power_table(x)
  }
  plot.data <- x
  
  
  # extract renamed power definition
  powerType <- NULL
  yLabel <- "power"
  
  
  # For M=1, we only have one kind of power.
  if ( M == 1 ) {
    stopifnot( is.null( power.definition ) || power.definition == "D1indiv" || power.definition == "indiv.mean" )
    power.definition = NULL
    yLabel = "individual power"
    powerType = "individual"
  } else if ( !is.null( power.definition ) ) {
    
    power.names <- get_power_names( M, long = TRUE )
    if ( !( power.definition %in% names(power.names) ) ) {
      sstop( "Power %s not one of %s",
             power.definition, paste0( names(power.names), collapse=", " ) )
    }
    powerType <- power.names[[power.definition]]

    pstat = parse_power_definition(power.definition)
    if( pstat$indiv ) {
      powerType <- "mean individual"
    }
    
    # filter to only relevant power definition
    plot.data <- plot.data %>%
      dplyr::filter(.data$power == powerType)
    
    yLabel <- paste0(powerType, " power")
  } else {
    if ( M > 1 ) {
      # drop individual power, just report the mean power.
      plot.data <- plot.data %>%
        dplyr::filter( !stringr::str_detect( .data$power, "individual outcome" ) )
    }
  }
  
  
  # pivot to long table, one row per MTP
  plot.data <- plot.data %>%
    dplyr::rename(powerType = .data$power) %>%
    tidyr::pivot_longer( cols = tidyselect::all_of( MTPs ),
                         names_to = "MTP", values_to = "power")
  
  # Aggregate data, if multiple things varying
  var_names = attr( x, "var_names" )
  if ( length( var_names ) > 1 ) {
      plot.data <- plot.data %>%
        dplyr::group_by( dplyr::across( c( "powerType", "MTP", var.vary ) ) ) %>%
        dplyr::summarise( power = mean( power ) )
      
      smessage('Note: Averaged power across other varying factors in grid: %s',
             paste0( setdiff( var_names, var.vary ), collapse = ", " ) )
  }
  
  
  # convert to factors for plotting
  plot.data <- plot.data %>%
    dplyr::mutate(target.power = as.numeric(.data$power),
                  MTP = as.factor(.data$MTP),
                  powerType = as.factor(.data$powerType))
  plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])
  
  
  # remove NA values
  plot.data <- plot.data %>%
    dplyr::filter(!is.na(.data$power))
 
  
  grid.plot <- ggplot2::ggplot(
    data = plot.data,
    ggplot2::aes_string(x = var.vary,
                        y = "power",
                        shape = "MTP",
                        color = "MTP")) +
    ggplot2::geom_point(size = 2,
                        position = ggplot2::position_dodge(width = 0.125)) +
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::ggtitle(paste0(powerType , " power when ", var.vary, " varies")) +
    ggplot2::labs(x = var.vary, y = yLabel,
                  color = "", shape = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                                      face = "bold",
                                                      vjust = 1,
                                                      hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 10))
  
  if ( is.null( power.definition) && M > 1 ) {
    grid.plot <- grid.plot + ggplot2::facet_wrap( ~ powerType )
  }
  
  return(grid.plot)
  
}





##### MDES and Sample grid plot methods #####

fetch_power_type = function( x, power.definition ) {
  M <- params(x)$M
  if ( !is.numeric(M) ) {
    stopifnot( !is.null( x$M ) )
    M = max( x$M )
  } 
  
  # extract renamed power definition
  power.names <- get_power_names( M, long = TRUE )
  powerType <- power.names[[power.definition]]
  
  return( powerType )
}


# # extract renamed power definition
# power.names <- get_power_names( M, long = TRUE )
# if (M != 1)
# {
#   powerType <- power.names[[power.definition]]
# } else
# {
#   powerType <- "individual outcome 1"
# }
# # rename for nicer graphing
# if(powerType == "individual power"){
#   powerType <- "individual"
# }




handle_power_definition = function( x, power.definition, outcome, var.vary ) {
  plot.data <- as.data.frame(x)
  
  multiPower = FALSE
  if (  "power.definition" %in% colnames(plot.data) )  {
    if ( is.null( power.definition ) ) {
      multiPower <- TRUE
    } else {
      plot.data = dplyr::filter( plot.data, .data$power.definition == power.definition )
      stopifnot( nrow(x) > 0 )
    }
  } else {
    if ( is.null( power.definition )  ) {
      power.definition <- params(x)$power.definition
    } else {
      stopifnot( (params(x)$power.definition == power.definition ) )
    }
  }
  
  powerType = ""
  if ( ! multiPower ) {
    powerType <- fetch_power_type( x, power.definition )
  } 
  
  powerPer = ""
  if ( is.numeric( params(x)$target.power ) ) {
    powerPer = paste0( round( 100 * params(x)$target.power ), "% " )
  } else {
    # target power is varying.  Do nothing and hope.
  }
  
  # converting data type for graphing purposes
  # plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])
  
  title <- ""
  if ( powerType == "" ) {
    title = paste0( outcome, " for ", powerPer, powerType,
                    " power when ", var.vary, " varies")
  } else {
    title = paste0( outcome, " for ", powerPer,
                    "power when ", var.vary, " varies")
  }
  
  list( plot.data = plot.data, powerType = powerType, multiPower = multiPower,
        title = title )
}




#' Plot a grid pump mdes object
#'
#' @inheritParams plot.pumpgridresult
#' @keywords internal
plot.pumpgridresult.mdes <- function( x, power.definition = NULL, var.vary, ...  )
{
  M <- params(x)$M
  
  res <- handle_power_definition( x, power.definition, outcome = "MDES", var.vary = var.vary )
  
  plot.data = res$plot.data
  
  # Aggregate data, if multiple things varying
  var_names = setdiff( attr( x, "var_names" ), c( "MTP", "power.definition" ) )
  if ( length( var_names ) > 1 ) {
    plot.data <- plot.data %>%
      dplyr::group_by( dplyr::across( c( "MTP", var.vary ) ) ) %>%
      dplyr::summarise( Adjusted.MDES = mean( Adjusted.MDES ) )
    
    smessage('Note: Averaged Adjusted.MDES across other varying factors in grid: %s',
             paste0( setdiff( var_names, var.vary ), collapse = ", " ) )
  }
  
  
  # converting data type for graphing purposes
  plot.data <- plot.data %>%
    dplyr::mutate(Adjusted.MDES = as.numeric(.data$Adjusted.MDES))
  plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])
  
  
  grid.plot <- ggplot2::ggplot(
    data = plot.data,
    ggplot2::aes_string(x = var.vary,
                        y = "Adjusted.MDES",
                        color = "MTP",
                        shape = "MTP")) +
    ggplot2::geom_point(size = 2,
                        position = ggplot2::position_dodge(width = 0.125)) +
    ggplot2::ggtitle( res$title ) + 
    ggplot2::labs(x = paste0(var.vary, " (same across all outcomes)"),
                  y = "MDES",
                  color = "",
                  shape = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                                      face = "bold",
                                                      vjust = 1,
                                                      hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 10))
  
  if ( res$multiPower ) {
    grid.plot <- grid.plot + ggplot2::facet_wrap( ~ power.definition )
  }
  
  return( grid.plot )
}



#' Plot a grid pump sample object
#'
#' @inheritParams plot.pumpgridresult
#' @keywords internal
plot.pumpgridresult.sample <- function( x, power.definition = NULL, var.vary, ...  ) {
  
  params = params(x)

  sampleType <- attr( x, "sample.level" )
  
  res <- handle_power_definition( x, power.definition, outcome = sampleType,
                                  var.vary = var.vary )
  
  plot.data = res$plot.data
  
  # Aggregate data, if multiple things varying
  var_names = setdiff( attr( x, "var_names" ), c( "MTP", "power.definition" ) )
  if ( length( var_names ) > 1 ) {
    plot.data <- plot.data %>%
      dplyr::group_by( dplyr::across( c( "MTP", var.vary ) ) ) %>%
      dplyr::summarise( Sample.size = mean( Sample.size ) )
    
    smessage('Note: Averaged Sample.size across other varying factors in grid: %s',
             paste0( setdiff( var_names, var.vary ), collapse = ", " ) )
  }
  
  # for nice axes
  if(max(plot.data$Sample.size) - min(plot.data$Sample.size) < 5)
  {
    ymin <- max(min(plot.data$Sample.size) - 3, 0)
    ymax <- max(plot.data$Sample.size) + 3
    
  } else
  {
    ymin <- min(plot.data$Sample.size)
    ymax <- max(plot.data$Sample.size)
  }
  integer.breaks <- function(ymin, ymax) {
      unique(floor(pretty(seq(ymin, (ymax + 1) * 1.1)))) 
  }
  
  
  grid.plot <- ggplot2::ggplot(
    data = plot.data,
    ggplot2::aes_string(x = var.vary,
                        y = "Sample.size",
                        color = "MTP",
                        shape = "MTP")) +
    ggplot2::geom_point(size = 2,
                        position = ggplot2::position_dodge(width = 0.125)) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax),
                                breaks = integer.breaks(ymin, ymax)) + 
    ggplot2::ggtitle( res$title ) +
    ggplot2::labs(x = paste0(var.vary, " (same across all outcomes)"),
                  y = "Sample size",
                  color = "",
                  shape = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                                      face = "bold",
                                                      vjust = 1,
                                                      hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_blank())
  
  if ( res$multiPower ) {
    grid.plot <- grid.plot + ggplot2::facet_wrap( ~ power.definition )
  }
  
  return( grid.plot )
}


#'@title Plot a pump grid result object
#'
#'@description
#'
#'Plots grid results across values of a single parameter, specified by the user
#'using var.vary, for a single definition of power, specified by
#'power.definition.
#'
#'If multiple things vary in the grid, the outcome (power, mdes, or sample size)
#'will be averaged (marginalized) across the other varying factors. This treats
#'the grid as a multifactor simulation, with this showing the "main effect" of
#'the specified parameter.
#'
#'@param x pumpgridresult object
#'@param power.definition definition of power to plot.  If NULL, plot all
#'  definitions as a facet wrap.
#'@param var.vary variable to vary on X axis.  If NULL, and only one thing
#'  varies, then it will default to single varying parameter.
#'@param ... additional parameters
#'
#'@return a ggplot object
#'@export
#'
#' @examples
#'g <- pump_power_grid( d_m = "d3.2_m3ff2rc", MTP = c( "HO", "BF" ),
#'  MDES = 0.10, J = seq(5, 10, 1), M = 5, K = 7, nbar = 58,
#'  Tbar = 0.50, alpha = 0.15, numCovar.1 = 1,
#'  numCovar.2 = 1, R2.1 = 0.1, R2.2 = 0.7,
#'  ICC.2 = 0.25, ICC.3 = 0.25, rho = 0.4, tnum = 500)
#'plot(g, power.definition = 'min1')

plot.pumpgridresult <- function( x, power.definition = NULL, var.vary = NULL, ... )
{
  # validation
  stopifnot( is.pumpgridresult( x ) )
  
  var_names <- setdiff( attr(x, "var_names" ), c( "MTP", "power.definition" ) )
  
  if ( is.null( var_names ) ) {
    stop( "No list of varying design elements found in pump grid result" )
  }
  
  if( !is.null( var.vary ) ) {
    if ( !(var.vary %in% var_names) ) {
      sstop('Please provide a var.vary amongst the variables that vary. %s is not listed.', var.vary )
    }
  } else {
    if ( length( var_names ) > 1 ) {
      sstop( "Need to specify what variable (of %s) to examine via var.vary parameter in plot()",
             paste0( var_names, collapse="," ) )
    }
    var.vary = var_names
  }
  
  
  if(pump_type(x) == 'power') {
    
    grid.plot <- plot.pumpgridresult.power(
        x, power.definition = power.definition,
        var.vary = var.vary, ... )
    
  } else if (pump_type(x) == 'mdes') {
    
    grid.plot <- plot.pumpgridresult.mdes(
        x, power.definition = power.definition,
        var.vary = var.vary, ... )
    
  } else if(pump_type(x) == 'sample') {
    
    grid.plot <- plot.pumpgridresult.sample(
        x, power.definition = power.definition,
        var.vary = var.vary, ... )
    
  } else {
    stop('Invalid pumpresult type.')
  }
  
  return(grid.plot)
}

