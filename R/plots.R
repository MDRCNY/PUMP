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
  if(is.na(pwr$Sample.size))
  {
    stop('plot_power_curve() does not work if optimizer has not converged. Try plot_power_search() for additional information.')
  }
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
  #lims <- grDevices::extendrange( r = range( test.pts$power, test.pts$target.power[[1]], na.rm = TRUE ), 0.15 )
  limsX <- grDevices::extendrange( r = xrng, 0.15 )

  plot1 <-  ggplot2::ggplot( test.pts ) +
    ggplot2::geom_hline( yintercept = test.pts$target.power[[1]], col = "purple" ) +
    ggplot2::theme_minimal() +
    ggplot2::stat_function( col="red", fun = function(x) { bounded_logistic_curve( x, params = fit ) } ) +
    ggplot2::guides(colour="none", size="none") +
    ggplot2::coord_cartesian( ylim=c(0,1), xlim = limsX ) +
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

#' Plot a single scenario pump object
#'
#' @param x pumpresult object
#'
#' @export
plot.pumpresult <- function( x )
{
  stopifnot( is.pumpresult( x ) )

  if(pump_type(x) == 'power')
  {
    if(attr( x, "long.table" ))
    {
      stop('Please call on a pump object not in long.table format')
    }

    plot.data <-
      x %>%
      dplyr::select_all() %>%
      dplyr::select(-indiv.mean) %>%
      tidyr::pivot_longer(!MTP, names_to = "powerType", values_to = "power")

    mtpname <- levels(as.factor(plot.data$MTP))[1]

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

    # single scenario plot
    ss.plot <- ggplot2::ggplot(
      data = plot.data,
      ggplot2::aes(x = powerType,
          y = power,
          shape = MTP,
          colour = MTP)) +
      ggplot2::geom_point(size = 1.5, position = ggplot2::position_dodge(0.25)) +
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
      ggplot2::labs(colour = "MTP", shape = "")
  } else if(power.type(x) == 'mdes')
  {
    stop('Not yet implemented')
  }  else if(power.type(x) == 'sample')
  {
    stop('Not yet implemented')
  } else
  {
    stop('Invalid pumpresult type.')
  }

  return(ss.plot)
}


#' Plot a pump grid object
#'
#' @param x pumpgrid object
#'
#' @export
plot.pumpgrid <- function( x, power.definition, var.vary )
{
  stopifnot( is.pumpgrid( x ) )

  if(pump_type(x) == 'power')
  {
    if(!attr( x, "long.table" ))
    {
      stop('Please call on a pump object in long.table format.')
    }

    if(ncol(x) > 6)
    {
      stop('Grid plot only works on grid objects with a single parameter varying.')
    }

    M <- params(x)$M
    res <- c('mean individual', 'complete',  paste(1:(M-1),'minimum', sep = '-'))
    names(res) <- c('indiv.mean', 'complete', paste('min', 1:(M-1), sep = ''))

    # Pulling out the power definition of interest matched with what's in the output table
    def_power_filter <- res[[power.definition]]

    # Pulling out only that power definition
    plot.data <-
      x %>%
      dplyr::filter(power %in% def_power_filter) %>%
      dplyr::mutate(power = ifelse(power == "mean individual", "individual power", power))

    plot.data <- plot.data %>%
      dplyr::select(-MDES) %>%
      dplyr::relocate(design)

    plot.data <-
      plot.data %>%
      dplyr::select_all() %>%
      dplyr::select(-design) %>%
      dplyr::arrange(desc(power)) %>%
      dplyr::rename(powerType = power) %>%
      tidyr::pivot_longer(!c(var.vary, powerType), names_to = "MTP", values_to = "power") %>%
      dplyr::filter(!stringr::str_detect(powerType, "individual outcome")) %>%
      dplyr::mutate(powerType = ifelse(MTP == "None",
                                       "raw mean individual",
                                       powerType))

    # converting Power Type to a factor for coloring
    plot.data$powerType <- factor(plot.data$powerType)

    # converting data type for graphing purposes
    plot.data <- plot.data %>%
      dplyr::mutate(target.power = as.numeric(power),
                    MTP = as.factor(MTP))

    # Converting to factor the variable that we are varying
    plot.data[[1]] <- as.factor(plot.data[[1]])

    # name of MTP
    mtpname <- levels(as.factor(plot.data$MTP))[1]
    powerType <- levels(as.factor(plot.data$powerType))[1]

    if(powerType == "individual power"){
      powerType <- "individual"
    }

    grid.plot <- ggplot2::ggplot(
      data = plot.data,
      ggplot2::aes_string(x = var.vary,
                 y = "power",
                 shape = "MTP",
                 color = "MTP")) +
      ggplot2::geom_point(size = 1.5, position = ggplot2::position_dodge(width = 0.125)) +
      ggplot2::scale_y_continuous(limits = c(0,1)) +
      ggplot2::ggtitle(paste0(powerType , " power when ", var.vary, " varies")) +
      ggplot2::labs(x = paste0(var.vary, " (constant across all outcomes)"),
           y = paste0(powerType, " power"),
           color = "",
           shape = "") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 16,
                                      face = "bold",
                                      vjust = 1,
                                      hjust = 0.5),
            axis.text = ggplot2::element_text(size = 10))

  } else if(power.type(x) == 'mdes')
  {
    stop('Not yet implemented')
  }  else if(power.type(x) == 'sample')
  {
    stop('Not yet implemented')
  } else
  {
    stop('Invalid pumpresult type.')
  }

  return(grid.plot)
}

