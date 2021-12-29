#' Examine search path of the power search.
#'
#' This will give a plot of power vs. mdes or sample size. It can be useful to
#' see how quickly power changes as a function of these design parameters. Can
#' also be useful to gauge where convergence went poorly.
#'
#' @param pwr Result from the pump_sample or pump_mdes (or data frame from,
#'   e.g., power_curve()).
#' @param plot.points flag; whether to plot individually tested points on curve
#' @param fit A four parameter bounded logstic curve (if NULL will fit one to passed
#'   points).
#'
#' @inheritParams power_curve
#'
#' @export
plot_power_curve <- function( pwr, plot.points = TRUE,
                              all = FALSE,
                              low = NULL, high = NULL, grid.size = 5, tnum = 2000,
                              fit = NULL ) {
  

  if ( is.pumpresult( pwr ) ) {
    if( is.na( pwr$Sample.size ) )
    {
        stop('plot_power_curve() does not work if algorithm has not converged.')
    }
    test.pts <- power_curve(pwr, low=low, high=high, grid.size = grid.size,
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
  
  xrng = range( test.pts$pt )
  #lims <- grDevices::extendrange( r = range( test.pts$power, test.pts$target.power[[1]], na.rm = TRUE ), 0.15 )
  limsX <- grDevices::extendrange( r = xrng, 0.15 )
  
  # sample size or MDES?  If dataframe passed, we don't know what we are
  # plotting so we turn to the scale of X to determine the lower bound.
  if ( limsX[1] <= 0 ) {
    if ( xrng[[2]] > 10 ) {
      limsX[1] = 1
    } else {
      limsX[1] = 0.0001
    }
  }
  
  plot1 <-  ggplot2::ggplot( test.pts ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[[1]], col = "purple" ) +
    ggplot2::theme_minimal() +
    ggplot2::stat_function( col = "red", fun = function(x) { bounded_logistic_curve( x, params = fit ) } ) +
    ggplot2::guides(colour = "none", size="none") +
    ggplot2::coord_cartesian( ylim = c(0,1), xlim = limsX ) +
    ggplot2::labs( x = x_label, y = "power" )

  delrange = diff( xrng )
  if ( delrange < 2 ) {
    # Use normal scale.
  } else if ( delrange >= 2 && delrange <= 10 ) {
    # Tick marks for each sample size.
    xpt = seq( floor( xrng[[1]] ), ceiling( xrng[[2]] ) )
    plot1 = plot1 +  ggplot2::scale_x_continuous( breaks = xpt )
  } else {
    plot1 = plot1 +  ggplot2::scale_x_log10()
  }


  if ( plot.points ) {
    plot1 <- plot1 + ggplot2::geom_point( ggplot2::aes( .data$pt, .data$power, size = .data$w ), 
                                          alpha = 0.5 )
  }

  return( plot1 )

}





#' Examine search path of the power search.
#'
#' This will give a triple-plots about how the search narrowed down into the final
#' estimate.  Can be useful to gauge where convergence went poorly.
#'
#' @param pwr Result from a pump_sample or pump_mdes call.
#' @param fit A fitted curve to the search.
#' 
#' @return a ggplot plot (a gridExtra arrangement of 3 plots, technically).
#'
#' @export
plot_power_search <- function( pwr, fit = NULL, target_line = NULL) {
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
  
  lims <- grDevices::extendrange( r = range( test.pts$power, test.pts$target.power[[1]], na.rm=TRUE ), 0.15 )
  
  plot1 <-  ggplot2::ggplot( test.pts ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[1], col = "purple" ) +
    ggplot2::geom_point( ggplot2::aes( .data$pt, .data$power, size = .data$w ), alpha = 0.5 ) +
    ggplot2::theme_minimal() +
    ggplot2::geom_text( ggplot2::aes( .data$pt, .data$power, label = .data$step ), vjust = "bottom", nudge_y = 0.01, size=3 ) +
    ggplot2::stat_function( col = "red", fun = function(x) { bounded_logistic_curve( x, params = fit ) } ) +
    ggplot2::guides(colour = "none", size="none") +
    ggplot2::coord_cartesian(ylim = lims ) +
    ggplot2::scale_x_log10()

  if ( !is.null( target_line ) ) {
    plot1 <- plot1 + ggplot2::geom_vline( xintercept = target_line, col = "purple" )
  }
  

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
  
  if ( !is.null( target_line ) ) {
    plot3 <- plot3 + ggplot2::geom_hline( yintercept = target_line, col = "purple" )
  }
    

  gridExtra::grid.arrange(plot1, plot2, plot3, ncol=3) #+
  #title( "Search path for optimize_power" )

}

#' Plot a single scenario pump object
#'
#' @param x pumpresult object
#' @param ... Additional parameters
#'
#' @export
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
      tidyr::pivot_longer(!.data$MTP, names_to = "powerType", values_to = "power")

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
      ggplot2::ggtitle(paste0("Adjusted power across different definitions of power")) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                      face = "bold",
                                      vjust = 1,
                                      hjust = 0.5),
            axis.text.x = ggplot2::element_text(size = 10, angle = 45),
            axis.text.y = ggplot2::element_text(size = 10),
            axis.title  = ggplot2::element_text(size = 10)
      ) +
      ggplot2::labs(color = "", shape = "")
  } else if( power.type(x) %in% c('mdes', 'sample') )
  {
    stop('plot() only works on pump_power() objects or grid objects, not pump_mdes() or pump_sample() objects.')
  } else
  {
    stop('Invalid pumpresult type.')
  }

  return(ss.plot)
}



#' Plot a grid pump power object
#'
#' @param x pumpgridresult object
#' @param power.definition user must choose one definition of power to plot
#' @param var.vary a single parameter to show varying on the plot
#' @param ... Additional parameters
#'
#' @export
plot.pumpgridresult.power <- function( x, power.definition, var.vary, ... ) {
  
  M <- params(x)$M
  MTPs <- unique(c("None", params(x)$MTP))
  
  if(!attr( x, "long.table" ))
  {
    x <- transpose_power_table(x)
  }
 
  # extract renamed power definition
  power.names <- get_power_names( M, long = TRUE )
  if (M != 1)
  {
    powerType <- power.names[[power.definition]]
  } else
  {
    powerType <- "individual outcome 1"
  }

  # filter to only relevant power definition
  plot.data <-
    x %>%
    dplyr::filter(.data$power == powerType)

  # remove MDES unless we want to vary it
  if(var.vary != 'MDES')
  {
    plot.data <- plot.data %>%
      dplyr::select(-.data$MDES)
  }
  
  # pivot to long table
  plot.data <-
      plot.data %>%
      dplyr::select_all() %>%
      dplyr::select(-.data$design) %>%
      dplyr::rename(powerType = .data$power) %>%
      tidyr::pivot_longer( cols = all_of( MTPs ),
                           names_to = "MTP", values_to = "power")

  # convert to factors for plotting
  plot.data <- plot.data %>%
      dplyr::mutate(target.power = as.numeric(.data$power),
                    MTP = as.factor(.data$MTP),
                    powerType = as.factor(.data$powerType))
  plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])
  
  # for pretty graph labels
  if(powerType == "individual power"){
    powerType <- "individual"
  }

  # remove NA values
  plot.data <- plot.data %>%
      dplyr::filter(!is.na(.data$power))
  
  # if two variables vary, send warning
  if(ncol(plot.data) > 4)
  {
      message('Note: more than one parameter varying in grid object.')
  }
  
  grid.plot <- ggplot2::ggplot(
    data = plot.data,
    ggplot2::aes_string(x = var.vary,
                        y = "power",
                        shape = "MTP",
                        color = "MTP")) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.125)) +
    ggplot2::scale_y_continuous(limits = c(0,1)) +
    ggplot2::ggtitle(paste0(powerType , " power when ", var.vary, " varies")) +
    ggplot2::labs(x = var.vary,
                  y = paste0(powerType, " power"),
                  color = "",
                  shape = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                                      face = "bold",
                                                      vjust = 1,
                                                      hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 10))
  
  return(grid.plot)
  
}


#' Plot a grid pump mdes object
#'
#' @inheritParams plot.pumpgridresult.power
#'
#' @export
plot.pumpgridresult.mdes <- function( x, power.definition, var.vary, ...  )
{
  M <- params(x)$M
  
  # extract renamed power definition
  power.names <- get_power_names( M, long = TRUE )
  if (M != 1)
  {
    powerType <- power.names[[power.definition]]
  } else
  {
    powerType <- "individual outcome 1"
  }
  # rename for nicer graphing
  if(powerType == "individual power"){
      powerType <- "individual"
  }
    
  # extract column name from table
  powerColName <- grep(power.definition, colnames(x), value = TRUE)

  # converting data type for graphing purposes
  plot.data <- x %>%
    dplyr::mutate(Adjusted.MDES = as.numeric(Adjusted.MDES))
  plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])

  # if two variables vary, send warning
  if(ncol(plot.data) > 6)
  {
      message('Note: more than one parameter varying in grid object.')
  }
  
  grid.plot <- ggplot2::ggplot(
    data = plot.data,
    ggplot2::aes_string(x = var.vary,
                        y = "Adjusted.MDES",
                        color = "MTP",
                        shape = "MTP")) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.125)) +
    ggplot2::ggtitle(paste0("MDES for ", powerType, " power when ", var.vary, " varies")) + 
    ggplot2::labs(x = paste0(var.vary, " (same across all outcomes)"),
                  y = "MDES",
                  color = "",
                  shape = "") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                                      face = "bold",
                                                      vjust = 1,
                                                      hjust = 0.5),
                   axis.text = ggplot2::element_text(size = 10))
  
  return( grid.plot )
}



#' Plot a grid pump sample object
#'
#' @inheritParams plot.pumpgridresult.power
#'
#' @export
plot.pumpgridresult.sample <- function( x, power.definition, var.vary, ...  ) {
  
  M <- params(x)$M
  sampleType <- attr( x, "sample.level" )

  # extract renamed power definition
  power.names <- get_power_names( M, long = TRUE )
  if (M != 1)
  {
    powerType <- power.names[[power.definition]]
  } else
  {
    powerType <- "individual outcome 1"
  }
  # rename for nicer graphing
  if(powerType == "individual power"){
      powerType <- "individual"
  }

  # extract column name from table
  powerColName <- grep(power.definition, colnames(x), value = TRUE)
  
  # converting data type for graphing purposes
  plot.data <- x %>%
    dplyr::mutate(Sample.size = as.numeric(Sample.size))
  plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])

  # if two variables vary, send warning
  if(ncol(plot.data) > 7)
  {
    message('Note: more than one parameter varying in grid object.')
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
  integer.breaks <- function(ymin, ymax) { unique(floor(pretty(seq(ymin, (ymax + 1) * 1.1)))) }
  
  grid.plot <- ggplot2::ggplot(
    data = plot.data,
    ggplot2::aes_string(x = var.vary,
                        y = "Sample.size",
                        color = "MTP",
                        shape = "MTP")) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.125)) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax),
                                breaks = integer.breaks(ymin, ymax)) + 
    ggplot2::ggtitle(paste0(sampleType, " for ", powerType, " power when ", var.vary, " varies")) + 
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
  
  return( grid.plot )
}


#' Plot a pump grid result object
#'
#' @param x pumpgridresult object
#' @param power.definition Definition of power to plot
#' @param var.vary Variable to vary on X axis
#'
#' @export
plot.pumpgridresult <- function( x, power.definition, var.vary, ... )
{
  # validation
  stopifnot( is.pumpgridresult( x ) )
  
  var_names <- attr(x, "var_names" )
  if ( is.null( var_names ) ) {
    stop( "No list of varying design elements found in pump grid result" )
  }
  if( !(var.vary %in% var_names) )
  {
      stop('Please provide a var.vary amongst the variables that vary')
  }

  if(pump_type(x) == 'power') {
  
    grid.plot <- plot.pumpgridresult.power(x, power.definition = power.definition,
                                           var.vary = var.vary, ... )
    
  } else if (pump_type(x) == 'mdes') {
     
    grid.plot <- plot.pumpgridresult.mdes(x, power.definition = power.definition,
                                          var.vary = var.vary, ... )
    
  } else if(pump_type(x) == 'sample') {
      
    grid.plot <- plot.pumpgridresult.sample(x, power.definition = power.definition,
                                            var.vary = var.vary, ... )
      
  } else {
    stop('Invalid pumpresult type.')
  }

  return(grid.plot)
}

