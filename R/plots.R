#' @title Examine a power curve (result function)
#'
#' @description This will give a plot of power vs. MDES or sample size. It can
#'   be useful to see how quickly power changes as a function of these design
#'   parameters. Can be useful to diagnose relatively flat power curves, where
#'   power changes little as a function of MDES or sample size, and can also be
#'   useful to gauge where convergence went poorly.
#'
#' @param pwr pumpresult object or data.frame; result from calling 
#' pump_sample or pump_mdes (or data frame from, e.g., power_curve()).
#' @param plot.points logical; whether to plot individually 
#' tested points on curve.
#' @param fit a four parameter bounded logistic curve 
#' (if NULL will fit one to passed points).
#' @param breaks scalar; the desired number of tick marks on the axes.
#' @inheritParams power_curve
#' 
#' @return plot; a ggplot object of power across values.
#' @keywords internal
plot_power_curve <- function(pwr, plot.points = TRUE,
                             all = TRUE,
                             low = NULL, high = NULL,
                             grid.size = 5, tnum = 2000,
                             breaks = grid.size, 
                             fit = NULL ) {
    
    sample_axes <- TRUE
    
    if ( is.pumpresult( pwr ) ) {
        
        test.pts <- power_curve(pwr, low = low, high = high, 
                                grid.size = grid.size,
                                tnum = params(pwr)$tnum, all = all )
        x_label <- pump_type(pwr)
        sample_axes <- x_label == "sample"
        
    } else {
        stopifnot( is.data.frame(pwr) )
        test.pts <- pwr
        x_label <- "parameter"
        sample_axes <- max( test.pts$pt >= 2, na.rm = TRUE )
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
            limsX[1] <- 0
        }
    }
    
    plot1 <-  ggplot2::ggplot( test.pts ) +
        ggplot2::geom_hline( yintercept = test.pts$target.power[[1]],
                             col = "purple" ) +
        default_theme()
    
    if ( plot.points ) {
        plot1 <- plot1 +
            ggplot2::geom_point(
                ggplot2::aes( .data$pt, .data$power, size = .data$w ), 
                alpha = 0.5 )
    }
    
    plot1 <- plot1 +
        ggplot2::stat_function( 
            col = "red",
            fun = function(x) { bounded_logistic_curve( x, params = fit ) } ) +
        ggplot2::guides(colour = "none", size = "none") +
        ggplot2::coord_cartesian( ylim = c(0,1), xlim = limsX ) +
        ggplot2::labs( x = x_label, y = "power" )
    
    if ( sample_axes ) {
        scale <- get_sample_size_scale( test.pts$pt, breaks = breaks,
                                        include.points = plot.points )
        plot1 <- plot1 + scale
    }
    
  
    return( plot1 )
    
    
}





#' @title Examine search path of a power search (result function)
#'
#' @description This will give triple-plots about 
#' how the search narrowed down into the
#' final estimate.  Can be useful to gauge where 
#' convergence went poorly.
#'
#' @param pwr pumpresult object; result from a 
#' pump_sample or pump_mdes call.
#' @param fit a fitted curve to the search.
#' @param target.line scalar; if non-NULL, add a reference line 
#' for the true power (if known, e.g., from a pump_power call). 
#'
#' @return plot; a ggplot object
#'  (a ggpubr arrangement of 3 plots, technically) of the
#'  search path.
#'
#' @keywords internal
plot_power_search <- function(pwr, fit = NULL, target.line = NULL) {
    if ( is.pumpresult(pwr) ) {
        test.pts <- search_path(pwr)
    } else if ( is.data.frame(pwr) ) {
        test.pts <- pwr
    } else {
        test.pts <- pwr$test.pts
    }
    
    if (is.null(test.pts))
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
            ggplot2::geom_vline( xintercept = target.line, col = "darkgrey" )
    }
    
    
    lims <- grDevices::extendrange( r = range( test.pts$power, 
                                               test.pts$target.power[[1]], 
                                               na.rm = TRUE ), 0.15 )
    plot2 <-  ggplot2::ggplot( 
        test.pts,
        ggplot2::aes(.data$step, .data$power, size = .data$w) ) +
        ggplot2::geom_hline( yintercept = test.pts$target.power[1],
                             col = "darkgrey" ) +
        ggplot2::geom_point( alpha = 0.5 ) +
        ggplot2::scale_x_continuous( breaks = 0:max(test.pts$step) ) +
        default_theme() +
        ggplot2::coord_cartesian(ylim = lims ) +
        ggplot2::guides(colour = "none", size = "none")
    
    plot3 <-  ggplot2::ggplot( 
        test.pts,
        ggplot2::aes(.data$step, .data$pt, size = .data$w) ) +
        ggplot2::geom_point( alpha = 0.5 ) +
        ggplot2::scale_x_continuous( breaks = 0:max(test.pts$step) ) +
        default_theme() +
        ggplot2::guides(colour = "none", size = "none") +
        ggplot2::scale_y_log10()
    
    if ( !is.null( target.line ) ) {
        plot3 <- plot3 +
            ggplot2::geom_hline( yintercept = target.line, col = "darkgrey" )
    }
    
    ggpubr::ggarrange(plot1, plot2, plot3, 
                      ncol = 3, common.legend = TRUE, legend = "bottom" )
}


#' @title Plot a pumpresult object (result function)
#'
#' @description For the object returned by pump_power(), visualizes
#'   different definitions of power across MTPs. For the object
#'   returned by pump_mdes() or pump_sample(), plot a power curve as a
#'   function of MDES or sample size, respectively.  This latter call
#'   will calculate power over a passed range from low to high to
#'   generate this curve.
#'
#'   Several of the passed parameters only apply to the mdes or sample
#'   versions, and are for controlling the grid search and plot.
#'
#'   For pump_power, will include standard errors of uncertainty on
#'   calculated power. These depend on number of iterations (tnum)
#'   used in the simulation.
#'
#' @param x pumpresult object.
#' @param type string; "power" or "search". Specifies whether to plot
#'   the default power graph, or the search path. The search path is
#'   only valid for MDES and SS results.
#' @param low Low range of x-axis and curve calculation for sample or
#'   MDES plots.  (Optional.)
#' @param high High range of x-axis and curve calculation. (Optional.)
#' @param all Logical. If TRUE, merge in the search path from the
#'   original search to the estimated power curve, for MDES or sample
#'   plots.
#' @param grid.size If calculating curve for sample or MDES plot, how
#'   many grid points?
#' @param breaks If plotting a curve for sample or MDES, where to put
#'   the grid points?
#' @param include_SE Include (approximate) SEs on the power estimates,
#'   if they are naturally calculated.
#' @param ... additional parameters, such as, in case of sample or
#'   mdes objects, tnum for setting number of replicates or all
#'   (logical) for determining whether to include original points in
#'   the estimated curve, or include.points  (logical) for including
#'   points on the plot itself.
#'
#' @return plot; a ggplot object of power across different
#'   definitions.
#'
#' @export
#'
#' @examples
#' pp1 <- pump_power(d_m = "d2.2_m2rc", MTP = 'HO',
#'  nbar = 50, J = 20, M = 8, numZero = 5,
#'  MDES = 0.30, Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
#'  numCovar.1 = 1, numCovar.2 = 1, R2.1 = 0.1, R2.2 = 0.7,
#'  ICC.2 = 0.05, rho = 0.2, tnum = 200)
#'
#' plot(pp1)
#'
#' J <- pump_sample(d_m = "d2.1_m2fc",
#'    MTP = 'HO', power.definition = 'D1indiv',
#'    typesample = 'J', target.power = 0.6,
#'    nbar = 50, M = 3, MDES = 0.125,
#'    Tbar = 0.5, alpha = 0.05,
#'    numCovar.1 = 1, R2.1 = 0.1, ICC.2 = 0.05,
#'    rho = 0.2, tnum = 500)
#' plot(J)
#' plot(J, type = "search")
#' 
plot.pumpresult <- function(x, type = "power", 
                            all = TRUE,
                            low = NULL, high = NULL,
                            grid.size = 5,
                            breaks = grid.size, 
                            include_SE = TRUE,
                            ... )
{
    stopifnot( is.pumpresult( x ) )
    stopifnot( type %in% c("power", "search") )
    
    if (pump_type(x) == "power")
    {
        if (type == "search" )
        {
            stop("Invalid plot type.")
        }
        if (attr( x, "long.table" ))
        {
            x <- transpose_power_table(x)
        }
        
        plot.data <-
            x %>%
            dplyr::select_all() %>%
            dplyr::select(-tidyselect::any_of("indiv.mean"), 
                          -tidyselect::starts_with("df"), 
                          -tidyselect::starts_with("SE") ) %>%
            tidyr::pivot_longer( !"MTP",
                                 names_to = "powerType", values_to = "power")
        
        # Creating power type as a factor for ordering on x axis
        M <- params(x)$M
        dPowers <- paste0("D", 1:M, "indiv")
        minPowers <- paste0("min", 1:M - 1)
        complete <- "complete"
        plot.data$powerType <- factor(
            plot.data$powerType,
            levels = c(dPowers, minPowers, complete),
            ordered = TRUE
        )
        
        # remove missing rows
        plot.data <- plot.data %>%
            as.data.frame( plot.data, row.names = NULL ) %>%
            dplyr::filter( !is.na( .data$power ) )
        
        plot.data <- dplyr::mutate( 
            plot.data, 
            SE = calc_binomial_SE( .data$power, params(x)$tnum ),
            CI_h = .data$power + 2 * .data$SE,
            CI_l = .data$power - 2 * .data$SE 
        )
        # single scenario plot
        ss.plot <- ggplot2::ggplot(
            data = plot.data,
            ggplot2::aes(x = .data$powerType,
                         y = .data$power,
                         shape = .data$MTP,
                         color = .data$MTP)) +
            ggplot2::geom_point( size = 2, 
                                 position = ggplot2::position_dodge(0.25) ) +
            ggplot2::scale_y_continuous(limits = c(0,1)) +
            ggplot2::ggtitle(paste0("Adjusted power across\n",
                                    "different definitions of power")) +
            default_theme() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, 
                                                               angle = 45),
                           axis.text.y = ggplot2::element_text(size = 10),
                           axis.title  = ggplot2::element_text(size = 10) ) +
            ggplot2::labs(color = "", shape = "")
        
        if ( include_SE ) {
            ss.plot = ss.plot + 
                ggplot2::geom_linerange( 
                    ggplot2::aes( ymin = .data$CI_l, ymax = .data$CI_h ),
                    position = ggplot2::position_dodge(0.25) )
        }
        return(ss.plot)
        
    } else if (pump_type(x) %in% c('mdes', 'sample') ) {
        
        if ( type == "power" ) {
            if ( is.null( low ) && pump_type( x ) == "mdes" ) {
                low <- 0  
            }
            return( plot_power_curve(x, low = low, high = high,
                                     grid.size = grid.size,
                                     breaks = breaks, ... ) )
        } else {
            return( plot_power_search(x, ... ) )
        }
    } else {
        stop('Invalid pumpresult type.')
    }
}



#' Plot a pump grid power object
#'
#' @inheritParams plot.pumpgridresult
#' 
#' @return Plot; a ggplot object
#' 
#' @importFrom stringr str_detect
#' @keywords internal
plot.pumpgridresult.power <- function(
        x, power.definition = NULL, var.vary = NULL, 
        color = "MTP",
        lines = TRUE, include.title = FALSE,  ... 
) {
    dots <- list( ... )
    
    M <- params(x)$M
    MTPs <- unique(c("None", params(x)$MTP))
    
    if (!attr( x, "long.table" )) {
        x <- transpose_power_table(x)
    }
    plot.data <- x
    
    # extract renamed power definition
    powerType <- NULL
    yLabel <- "power"
    
    # For M=1, we only have one kind of power.
    if ( M == 1 ) {
        stopifnot( is.null( power.definition ) || 
                       power.definition == "D1indiv" || 
                       power.definition == "indiv.mean" )
        power.definition <- NULL
        yLabel <- "individual power"
        powerType <- "individual"
    } else if ( !is.null( power.definition ) ) {
        
        power.names <- get_power_names( M, long = TRUE )
        if ( !( power.definition %in% names(power.names) ) ) {
            sstop( "Power %s not one of %s",
                   power.definition, 
                   paste0( names(power.names), collapse = ", " ) )
        }
        powerType <- power.names[[power.definition]]
        
        pstat <- parse_power_definition(power.definition)
        if ( pstat$indiv ) {
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
                dplyr::filter( 
                    !stringr::str_detect( .data$power, 
                                          "individual outcome" ) 
                )
        }
    }
    
    
    # pivot to long table, one row per MTP
    plot.data <- plot.data %>%
        dplyr::rename(powerType = "power") %>%
        tidyr::pivot_longer( cols = tidyselect::all_of( MTPs ),
                             names_to = "MTP", values_to = "power")
    
    # Aggregate data, if multiple things varying
    var_names <- setdiff( attr( x, "var_names" ), 
                          c( color, "power.definition" ) )
    
    if ( "MTP" %in% var_names ) {
        stop( "Cannot have multiple MTP for plotting if color not set to MTP" )
    }
    
    
    if ( length( var_names ) > 1 ) {
        plot.data <- plot.data %>%
            dplyr::group_by( 
                dplyr::across( tidyselect::all_of( 
                    c( "powerType", "MTP", var.vary ) ) 
                )
            ) %>%
            dplyr::summarise( power = mean( .data$power ) )
        
        smessage('Note: Averaged power across other varying 
                 factors in grid: %s',
                 paste0( setdiff( var_names, var.vary ), 
                         collapse = ", " ) )
    }
    
    
    # convert to factors for plotting

    # ensure our color is a factor
    plot.data$color <- as.factor( plot.data[[color]] )
    
    plot.data <- plot.data %>%
        dplyr::mutate(target.power = as.numeric(.data$power),
                      MTP = as.factor(.data$MTP),
                      powerType = as.factor(.data$powerType))
    if ( !is.numeric( plot.data[[var.vary]] ) ) {
        plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])
    } 
    
    # remove NA values
    plot.data <- plot.data %>%
        dplyr::filter(!is.na(.data$power))
    
    grid.plot <- ggplot2::ggplot(
        data = plot.data,
        ggplot2::aes(x = .data[[var.vary]],
                     y = .data$power,
                     shape = .data$color,
                     color = .data$color))
    
    if ( lines ) {
        grid.plot <- grid.plot +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_line()
    } else {
        grid.plot <- grid.plot + ggplot2::geom_point(
            size = 2,
            position = ggplot2::position_dodge(width = 0.125))
    }
    
    title <- NULL
    if ( include.title ) {
        title <- paste0(powerType , " power when ", var.vary, " varies")
    } 
    
    grid.plot <- grid.plot +
        ggplot2::scale_y_continuous(limits = c(0,1)) +
        ggplot2::labs(x = var.vary, y = yLabel, title = title,
                      color = color, shape = color) +
        default_theme()
    
    if ( lines ) {
        grid.plot <- grid.plot + ggplot2::geom_line()
    }
    
    if ( is.null( power.definition) && M > 1 ) {
        grid.plot <- grid.plot + 
            ggplot2::facet_wrap( ~ powerType, 
                                 nrow = dots$nrow, ncol = dots$ncol )
    }
    
    return(grid.plot)
    
}





##### MDES and Sample grid plot methods #####

fetch_power_type <- function(x, power.definition) {
    M <- params(x)$M
    if ( !is.numeric(M) ) {
        stopifnot( !is.null( x$M ) )
        M <- max( x$M )
    } 
    
    # extract renamed power definition
    power.names <- get_power_names( M, long = TRUE )
    powerType <- power.names[[power.definition]]
    
    return( powerType )
}


#' Plot a grid pump mdes object
#'
#' @inheritParams plot.pumpgridresult
#' @keywords internal
plot.pumpgridresult.mdes <- function( 
        x, power.definition = NULL, var.vary, 
        color = "MTP",
        lines = TRUE, include.title = FALSE, ...  
) {
    M <- params(x)$M
    dots <- list( ... )
    
    res <- handle_power_definition( 
        x, power.definition, outcome = "MDES", var.vary = var.vary, 
        include.title = include.title )
    
    plot.data <- res$plot.data
    
    # Aggregate data, if multiple things varying
    var_names <- setdiff( attr( x, "var_names" ), 
                          c( color, "power.definition" ) )
    
    if ( "MTP" %in% var_names ) {
        stop( "Cannot have multiple MTP for plotting if color not set to MTP" )
    }
    
    if ( length( var_names ) > 1 ) {
        plot.data <- plot.data %>%
            dplyr::group_by( dplyr::across( tidyselect::all_of( 
                c( "MTP", color, var.vary ) ) ) 
            ) %>%
            dplyr::summarise( Adjusted.MDES = mean( .data$Adjusted.MDES ) )
        
        smessage(paste('Note: Averaged Adjusted.MDES across other varying',
                       'factors in grid: %s'),
                 paste0( setdiff( var_names, var.vary ), collapse = ", " ) )
    }
    
    
    # converting data type for graphing purposes
    plot.data <- plot.data %>%
        dplyr::mutate(Adjusted.MDES = as.numeric(.data$Adjusted.MDES))
    
    # ensure our color is a factor
    plot.data$color <- as.factor( plot.data[[color]] )
    
    grid.plot <- ggplot2::ggplot(
        data = plot.data,
        ggplot2::aes(x = .data[[var.vary]],
                     y = .data$Adjusted.MDES,
                     color = .data$color,
                     shape = .data$color))
    
    if ( lines ) {
        grid.plot <- grid.plot +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_line()
    } else {
        grid.plot <- grid.plot + 
            ggplot2::geom_point(size = 2,
                                position = ggplot2::position_dodge(width = 0.125))
    }
    
    c_title <- color
    
    x_lab <- var.vary
    if ( M > 1 ) {
        x_lab = paste0(var.vary, " (same across all outcomes)")
    }
    
    grid.plot <- grid.plot +
        ggplot2::ggtitle( res$title ) + 
        ggplot2::labs(x = x_lab,
                      y = "MDES",
                      color = c_title,
                      shape = c_title) +
        default_theme() +
        ggplot2::expand_limits(y = 0)
    
    if ( lines ) {
        grid.plot <- grid.plot + ggplot2::geom_line()
    }
    
    if ( res$multiPower ) {
        grid.plot <- grid.plot + 
            ggplot2::facet_wrap( ~ power.definition, 
                                 nrow = dots$nrow, ncol = dots$ncol )
    }
    
    return( grid.plot )
}



#' Plot a grid pump sample object
#'
#' @inheritParams plot.pumpgridresult
#' @keywords internal
plot.pumpgridresult.sample <- function( 
        x, power.definition = NULL, var.vary, 
        color = "MTP",
        lines = TRUE, include.title = FALSE, ...  
) {
    dots <- list( ... )
    
    params <- params(x)
    
    sampleType <- attr( x, "sample.level" )
    
    res <- handle_power_definition( 
        x, power.definition, outcome = sampleType,
        include.title = include.title, var.vary = var.vary )
    
    plot.data <- res$plot.data
    
    # Aggregate data, if multiple things varying
    var_names <- setdiff( attr( x, "var_names" ), 
                          c( color, "power.definition" ) )
    if ( length( var_names ) > 1 ) {
        plot.data <- plot.data %>%
            dplyr::group_by( 
                dplyr::across( tidyselect::all_of( c( color, var.vary ) ) ) 
            ) %>%
            dplyr::summarise( Sample.size = 
                              mean( .data$Sample.size, na.rm = TRUE ) )
        
        smessage(paste('Note: Averaged Sample.size across other',
                       'varying factors in grid: %s'),
                 paste0( setdiff( var_names, var.vary ), collapse = ", " ) )
    }
    
    plot.data = dplyr::filter( plot.data, !is.na( .data$Sample.size ) )
    
    # for nice axes
    if (max(plot.data$Sample.size) - min(plot.data$Sample.size) < 5)
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
    
    # ensure our color is a factor
    plot.data$color <- as.factor( plot.data[[color]] )
    
    grid.plot <- ggplot2::ggplot(
        data = plot.data,
        ggplot2::aes(x = .data[[var.vary]],
                     y = .data$Sample.size,
                     color = .data$color,
                     shape = .data$color))
    
    if ( lines ) {
        grid.plot <- grid.plot +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_line()
    } else {
        grid.plot <- grid.plot + ggplot2::geom_point(
            size = 2,
            position = ggplot2::position_dodge(width = 0.125))
    }
    
    x_lab <- var.vary
    if ( params(x)$M > 1 ) {
        x_lab <- paste0(var.vary, " (same across all outcomes)")
    }
    
    grid.plot <- grid.plot + 
        ggplot2::scale_y_continuous(limits = c(ymin, ymax),
                                    breaks = integer.breaks(ymin, ymax)) + 
        ggplot2::ggtitle( res$title ) +
        ggplot2::labs(x = x_lab,
                      y = "Sample size",
                      color = "",
                      shape = "") +
        default_theme() +
        ggplot2::theme(legend.title = ggplot2::element_blank()) +
        ggplot2::expand_limits(y = 0)
    
    if ( res$multiPower ) {
        grid.plot <- grid.plot + 
            ggplot2::facet_wrap( ~ power.definition, 
                                 nrow = dots$nrow, ncol = dots$ncol )
    }
    
    return( grid.plot )
}


#'@title Plot a pumpgridresult object (result function)
#'
#'@description
#'
#'Plots grid results across values of a single parameter, specified by
#'the user using var.vary, for a single definition of power, specified
#'by power.definition.
#'
#'If multiple things vary in the grid, the outcome (power, mdes, or
#'sample size) will be averaged (marginalized) across the other
#'varying factors. This treats the grid as a multifactor simulation,
#'with this showing the "main effect" of the specified parameter.
#'
#'
#'@param x pumpgridresult object.
#'@param power.definition string; definition of power to plot. If
#'  NULL, plot all definitions as a facet wrap.
#'@param var.vary string; variable to vary on X axis. If NULL, and
#'  only one thing varies, then it will default to single varying
#'  parameter.
#'@param color string; Group lines by this element to make an
#'  interaction plot (default "MTP", giving one curve for each MTP).
#'@param lines logical; TRUE means connect dots with lines on the
#'  plots. FALSE means no lines.
#'@param include.title logical; whether to include/exclude title (if
#'  planning a facet wrap, for example).
#'@param ... additional parameters.
#'
#'@return plot; a ggplot object of outcome across parameter values.
#'
#'@export
#'
#' @examples
#'g <- pump_power_grid( d_m = "d3.2_m3ff2rc", MTP = c( "HO", "BF" ),
#'  MDES = 0.10, J = seq(5, 10, 1), M = 5, K = 7, nbar = 58,
#'  Tbar = 0.50, alpha = 0.15, numCovar.1 = 1,
#'  numCovar.2 = 1, R2.1 = 0.1, R2.2 = 0.7,
#'  ICC.2 = 0.25, ICC.3 = 0.25, rho = 0.4, tnum = 200)
#'plot(g, power.definition = 'min1')

plot.pumpgridresult <- function( 
        x, power.definition = NULL, var.vary = NULL, 
        color = "MTP",
        lines = TRUE, include.title = FALSE, ... 
)
{
    # validation
    stopifnot( is.pumpgridresult( x ) )
    
    var_names <- setdiff( attr(x, "var_names" ),
                          c( "MTP", color, "power.definition" ) )
    
    if ( is.null( var_names ) ) {
        stop( "No list of varying design elements found in pump grid result" )
    }
    
    if ( !is.null( var.vary ) ) {
        if ( !(var.vary %in% var_names) ) {
            sstop('Please provide a var.vary amongst the variables that vary. 
                  "%s" is not listed.', 
                  var.vary )
        }
    } else {
        if ( length( var_names ) > 1 ) {
            # Make a separate plot for each varying element!
            mps <- purrr::map( 
                var_names, plot.pumpgridresult, x = x, 
                power.definition = power.definition, color = color, 
                lines = lines, include.title = FALSE, ... )
            
            gd <- ggpubr::ggarrange( 
                plotlist = mps, common.legend = TRUE, 
                legend = "bottom", ncol = length(mps) ) 
            
            gd <- ggpubr::annotate_figure( 
                gd, top = paste0( "Main effects of varying ", 
                                  paste0( var_names, collapse = ", " ) ) )
            
            return( gd )
        } else {
            var.vary <- var_names
        }
    }
    
    
    if (pump_type(x) == 'power') {
        
        grid.plot <- plot.pumpgridresult.power(
            x, power.definition = power.definition, color = color,
            var.vary = var.vary, lines = lines, 
            include.title = include.title, ... )
        
    } else if (pump_type(x) == 'mdes') {
        
        grid.plot <- plot.pumpgridresult.mdes(
            x, power.definition = power.definition, color = color,
            var.vary = var.vary, lines = lines, 
            include.title = include.title, ... )
        
    } else if (pump_type(x) == 'sample') {
        
        grid.plot <- plot.pumpgridresult.sample(
            x, power.definition = power.definition, color = color,
            var.vary = var.vary, lines = lines, 
            include.title = include.title, ... )
        
    } else {
        stop('Invalid pumpresult type.')
    }
    
    return(grid.plot)
}





##### Helper functions for plotting #####

get_sample_tick_marks <- function( 
        desired_pts, breaks = 5, include.points = TRUE, log = FALSE 
) {
    
    desired_pts <- sort( unique( round( desired_pts ) ) )
    
    if ( log ) {
        desired_pts <- log( desired_pts )
    }
    
    mn <- min( desired_pts, na.rm = TRUE )
    mx <- max( desired_pts, na.rm = TRUE )
    rng <- mx - mn
    
    if ( log == FALSE && rng < breaks ) {
        breaks <- round(rng)
    }
    
    pts <- seq( mn, mx, length.out = breaks )
    
    if ( include.points ) {
        
        bw <- rng / (breaks + 1)
        
        rg <- cut( desired_pts, 
                   breaks = seq( mn - bw/2, mx + bw/2, 
                                 length.out = breaks + 1 ), 
                   labels = pts )
        
        grab_pt <- function(pt, gpts) {
            if ( length( gpts) == 0 ) {
                pt
            } else {
                dels <- abs( gpts - pt )
                gpts[ which.min( dels ) ]
            }
        }
        
        gps <- split( desired_pts, rg )
        
        vals <- purrr::map2_dbl( pts, gps, grab_pt )
        
        pts <- vals
    } 
    
    if ( log ) {
        return( unique( round( exp( pts ) ) ) )
    } else {
        return( round( pts ) )
    }
}






get_sample_size_scale <- function( 
        points, breaks = 5, include.points = FALSE 
) {
    
    delrange <- diff( range( points, na.rm = TRUE ) )
    if ( delrange > 50 ) {
        xpt <- get_sample_tick_marks( desired_pts = points, breaks = breaks, 
                                      include.points = include.points,
                                      log = TRUE )
        ggplot2::scale_x_log10( breaks = xpt )
    } else if ( delrange >= 2 && delrange <= 15 ) {
        # Tick marks for each sample size.
        xpt <- seq( round( min( points, na.rm = TRUE ) - 0.25 ), 
                    round( max( points, na.rm = TRUE ) + 0.25 ) )
        ggplot2::scale_x_continuous( breaks = xpt )
    } else {
        xpt <- get_sample_tick_marks( desired_pts = points, breaks = breaks, 
                                      include.points = include.points,
                                      log = FALSE )
        ggplot2::scale_x_continuous( breaks = xpt )
    }
}





handle_power_definition <- function( 
        x, power.definition, outcome, var.vary, include.title 
) {
    plot.data <- as.data.frame(x)
    
    multiPower <- FALSE
    if (  "power.definition" %in% colnames(plot.data) )  {
        if ( is.null( power.definition ) ) {
            multiPower <- TRUE
        } else {
            plot.data <- dplyr::filter( 
                plot.data, .data$power.definition == power.definition 
            )
            stopifnot( nrow(x) > 0 )
        }
    } else {
        if ( is.null( power.definition )  ) {
            power.definition <- params(x)$power.definition
        } else {
            stopifnot( (params(x)$power.definition == power.definition ) )
        }
    }
    
    powerType <- ""
    if ( !multiPower ) {
        powerType <- fetch_power_type( x, power.definition )
    } 
    
    powerPer <- ""
    if ( is.numeric( params(x)$target.power ) ) {
        powerPer <- paste0( round( 100 * params(x)$target.power ), "% " )
    } else {
        # target power is varying.  Do nothing and hope.
    }
    
    # converting data type for graphing purposes
    if ( !is.numeric( plot.data[[var.vary]] ) ) {
        plot.data[[var.vary]] <- as.factor(plot.data[[var.vary]])
    }
    
    title <- NULL
    if ( include.title ) {
        if ( powerType == "" ) {
            title <- paste0( outcome, " for ", powerPer, powerType,
                             " power when ", var.vary, " varies")
        } else {
            title <- paste0( outcome, " for ", powerPer,
                             "power when ", var.vary, " varies")
        }
    }
    
    list( plot.data = plot.data, powerType = powerType, 
          multiPower = multiPower,
          title = title )
}




default_theme <- function() {
    list( ggplot2::theme_minimal(),
          ggplot2::theme( legend.position = "bottom", 
                          legend.direction = "horizontal", 
                          legend.key.width = grid::unit(1,"cm"),
                          panel.border = ggplot2::element_blank(),
                          plot.title = ggplot2::element_text(size = 16,
                                                             face = "bold",
                                                             vjust = 1,
                                                             hjust = 0.5),
                          axis.text = ggplot2::element_text(size = 10)),
          ggthemes::scale_color_colorblind()
    )
}
