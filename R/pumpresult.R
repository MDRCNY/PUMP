
# exact: TRUE if direct calculation, FALSE if stocastic, with
# uncertainty.
make.pumpresult <- function(x,
                            type = c( "power", "mdes", "sample" ),
                            d_m = d_m,
                            params.list = NULL,
                            tries = NULL,
                            flat = FALSE,
                            exact = FALSE,
                            ...)
{
    type <- match.arg(type)
    class(x) <- c( "pumpresult", class(x) )
    attr(x, "type" ) <- type
    attr(x, "params.list") <- params.list
    attr(x, "d_m") <- d_m
    attr(x, "exact" ) <- exact
    ll <- list(...)
    for (l in names(ll)) {
        attr(x, l) <- ll[[ l ]]
    }
    if ( !is.null( tries ) ) {
        attr( x, "tries" ) <- tries
        attr( x, "search.range" ) <- c( min = min( tries$pt, na.rm = TRUE ),
                                        final = tries$pt[ nrow(tries) ],
                                        max = max( tries$pt, na.rm = TRUE ) )
        attr( x, "flat" ) <- flat
    }
    return( x )
}



#' @title Update a single pump call to a grid call (grid function)
#'
#' @description Take a pumpresult and provide lists 
#' of parameters to explore various versions 
#' of the initial scenario.
#' 
#'
#' @param x pump result object.
#' @param ... list of parameters to expand into a grid.
#' 
#' @return a pumpgridresult object; 
#' result of calling corresponding grid.
#'
#' @export
#' 
#' @examples
#' pp <- pump_power(d_m = "d2.1_m2fc", MTP = "HO",
#'   nbar = 200, J = 20, MDES = 0.2, M = 3,
#'   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
#'   R2.1 = 0.1, ICC.2 = 0.05, rho = 0, tnum = 500)
#'
#' gd <- update_grid( pp, J = c( 10, 20, 30 ) )
#'
update_grid <- function(x, ...)
{
    params <- attr(x,"param")
    params["d_m"] <- d_m(x)
    for (p in names(params)) {
        params[[p]] <- unique( params[[p]] )
    }
    pparam <- attr( x, "power.params.list" )
    params <- c( params, pparam )
    
    dts <- list(...)
    for (d in names(dts)) {
        params[[d]] <- dts[[d]]
    }
    result_type <- attr( x, "type" )
    if ( result_type == "power" ) {
        params$MDES <- unique(params[["MDES"]])
        do.call(pump_power_grid, params)
    } else if ( result_type == "mdes" ) {
        do.call( pump_mdes_grid, params )
    } else if ( result_type == "sample" ) {
        params["MDES"] <- unique(params[["MDES"]])
        if ( is.null( params[["typesample"]] ) ) {
            params["typesample"] <- attr( x, "sample.level" )
        }
        params[params$typesample] <- NULL
        do.call( pump_sample_grid, params )
    } else {
        stop(sprintf( "Unrecognized type, %s, in update_grid()", result_type ))
    }
}






#' @title Update a pump call, tweaking some parameters (core function)
#'
#' @description Works on objects returned by
#' pump_power(), pump_mdes(), or pump_sample().
#' One of the optional parameters can 
#' be a `type = something` argument, where
#' the "something" is either "power", "sample", or "mdes", 
#' if the call should be
#' shifted to a different pump call 
#' (pump_power, pump_sample, or pump_mdes,
#' respectively).
#'
#' @param object pump result object.
#' @param ... parameters as specified in `pump_power`, `pump_mdes`, and
#'   `pump_sample` that should be overwritten.
#' @param type string; can be "power", "mdes" or "sample", sets the 
#' type of the updated call (can be different from original).
#'
#' @return a pumpresult object: results of a new call 
#' using parameters of old object with newly 
#' specified parameters replaced.
#'
#' @export
#' 
#' @examples
#' ss <- pump_sample( d_m = "d2.1_m2fc", MTP = "HO",
#'   typesample = "J", nbar = 200, power.definition = "min1",
#'   M = 5, MDES = 0.05, target.power = 0.5, tol = 0.05,
#'   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, R2.1 = 0.1,
#'   ICC.2 = 0.15, rho = 0, final.tnum = 1000 )
#'
#' up <- update(ss, nbar = 40, tnum = 2000 )
#'
update.pumpresult <- function(object, type = NULL, ...)
{
    params <- params(object)
    orig_result_type <- attr(object, "type" )
    params["type"] <- orig_result_type
    params["d_m"] <- d_m(object)
    
    # for any vectors, collapse to single value of vector is identical
    # (this allows change in number of outcomes more easily)
    for (i in seq_along(params)) {
        if ( length( params[[i]] ) > 1 && 
             length( unique( params[[i]] ) ) == 1 )
        {
            params[i] = params[[i]][[1]]
        }
    }
    
    # Get new parameters
    dts <- list(...)
    
    # Are we changing what kind of calculation we want to perform?  
    # If so, adjust
    # some parameters as needed.
    
    # orig_result_type - the old type
    # type - the new type
    if ( !is.null(type) && type != orig_result_type ) {
        
        # Copy over sample size from the pump_sample call
        if ( orig_result_type == "sample" ) {
            ss <- object$`Sample.size`
            slvl <- attr(object, "sample.level" )
            params[[slvl]] <- ss
        }
        
        if ( type == "power" ) {
            params["target.power"] <- NULL
            params["power.definition"] <- NULL
            params["tol"] <- NULL
            params["start.tnum"] <- NULL
            params["final.tnum"] <- NULL
            params["max.steps"] <- NULL
        }
        
        if ( type == "sample" ) {
            if ( !is.null( dts$typesample) ) {
                params[dts$typesample] <- NULL
            } else {
                stop( "Need to specify typesample for update to sample call" )
            }
        }
        
        if ( type == "mdes" ) {
            params["MDES"] <- NULL
        }
        
    } else {
        if ( orig_result_type == "sample" ) {
            params["typesample"] <- attr( object, "sample.level" )
        }
        type <- orig_result_type
    }
    
    for (d in names(dts)) {
        params[[d]] <- dts[[d]]
    }
    params$type <- NULL
    
    if ( type == "power" ) {
        do.call(pump_power, params)
    } else if ( type == "mdes" ) {
        do.call( pump_mdes, params )
    } else if ( type == "sample" ) {
        do.call( pump_sample, params )
    } else {
        stop( sprintf( "Unrecognized type, %s, in update()", type ) )
    }
}






#' @title pumpresult object for results of power calculations
#' @name pumpresult
#'
#' @description
#' The pumpresult object is an S3 class that holds the results from
#' `pump_power()`, 
#' `pump_sample()`, and `pump_mdes()`.
#'
#' It has several methods that pull different information from this object, and
#' some printing methods for getting nicely formatted results.
#'
#' Pump result objects are also data.frames, so they can be easily manipulated
#' and combined.  The return values from the `grid` functions will just return
#' data frames in general.
#'
#' @seealso update
#' @seealso update_grid
#' @seealso print_context
#' 
#' @param x a pumpresult object (except for is.pumpresult, 
#' where it is a generic object to check).
#' @rdname pumpresult
#' @examples
#' pp <- pump_power(d_m = "d3.2_m3ff2rc",
#'   MTP = 'HO', nbar = 50, J = 30, K = 10,
#'   M = 5, MDES = 0.125, Tbar = 0.5, alpha = 0.05,
#'   numCovar.1 = 1, numCovar.2 = 1,
#'   R2.1 = 0.1, R2.2 = 0.1, ICC.2 = 0.2, ICC.3 = 0.2,
#'   omega.2 = 0, omega.3 = 0.1, rho = 0.5, tnum = 1000)
#'   
#' print(pp)
#' params(pp)
#' print_context(pp)
#' d_m(pp)
#' pump_type(pp)
#' is.pumpresult(pp)
#' as.data.frame(pp)
#' dim(pp)
#' summary(pp)
#' transpose_power_table(pp)
#' 
#' J <- pump_sample(d_m = "d2.1_m2fc",
#'   MTP = 'HO', power.definition = 'D1indiv',
#'   typesample = 'J', target.power = 0.7,
#'   nbar = 50, M = 3, MDES = 0.125,
#'   Tbar = 0.5, alpha = 0.05, numCovar.1 = 1,
#'   R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2, tnum = 1000)
#'   
#' search_path(J)
#' power_curve(J)   
NULL




#' @title Get design and model parameters from 
#' pumpresult object (result function)
#'
#' @return params: List of design parameters used.
#'
#' @rdname pumpresult
#' @export
params <- function(x, ...)
{
    stopifnot( is.pumpresult( x ) || is.pumpgridresult( x ) )
    
    pp <- attr( x, "params.list" )
    pp_pow <- attr(x, "power.params.list" )
    if ( !is.null( pp_pow ) ) {
        pp <- c( pp, pp_pow )
    }
    return( pp )
}




#' @title Get context (design and model) from 
#' pumpresult object (result function)
#'
#' @return d_m: Context (d_m) used (as string).
#'
#' @rdname pumpresult
#' @export
#' 
d_m <- function(x, ...)
{
    stopifnot( is.pumpresult( x ) || is.pumpgridresult(x) )
    
    pp <- attr( x, "d_m" )
    return( pp )
}



#' @title Get design from pumpresult object (result function)
#'
#' @return design (the randomization and levels) as string.
#'
#' @rdname pumpresult
#' @export
#' 
design <- function(x, ...)
{
    pp = d_m( x )
    design = parse_d_m( pp )
    return( design$design )
}



#' @title Obtain full search path of pump_mdes or pump_sample call
#'
#' @return search_path: Dataframe describing search path, 
#' if it was saved in the pumpresult object.
#' @rdname pumpresult
#'
#' @export
#' 
search_path <- function(x, ...)
{
    stopifnot( is.pumpresult( x ) )
    rs <- attr( x, "tries" )
    if ( !is.null( rs ) ) {
        rs$delta <- rs$power - rs$target.power
    }
    return( rs )
}


# Was the calculation exact?
# 
# 
exact_calc <- function(x, ...)
{
    stopifnot( is.pumpresult( x ) )
    exact <- attr( x, "exact" )
    if ( is.null( exact ) ) {
        return( FALSE )
    } else {
        return( exact )
    }
}

#' @title Obtain a power curve for a range of sample size or MDES values
#'
#' @description This is used to see how power changes as a function of
#'   sample size or MDES.  It takes a fit pumpresult and calculates a power
#'   curve based on that scenario coupled with a passed range of
#'   values to make the curve over.
#'
#' @param x a pumpresult object.
#' @param all logical; if TRUE, merge in the search path from the
#'   original search.
#' @param low scalar; low range for curve.
#' @param high scalar; high range for the curve.
#' @param grid.size scalar; number of points to calculate power for.
#' @param tnum scalar; number of iterations to calculate power at each
#'   grid point.
#'
#' @importFrom rlang .data
#'
#' @return data.frame of power results.
#'
#' @export
#' 
power_curve <- function(x, all = FALSE,
                        low = NULL, high = NULL, 
                        grid.size = 5, tnum = 2000)
{
    stopifnot( is.pumpresult( x ) )
    fin_pts <- attr( x, "final.pts" )
    if ( is.null( fin_pts ) ) {
        fin_pts <- estimate_power_curve( x, 
                                         low = low, high = high, 
                                         grid.size = grid.size,
                                         tnum = tnum )
    }
    srch <- search_path( x )
    if ( all && !is.null( srch ) ) {
        srch <- dplyr::filter( srch, .data$pt < max( fin_pts$pt * 1.1 ) )
        fin_pts <- dplyr::bind_rows( fin_pts, srch ) %>%
            dplyr::arrange( .data$pt ) %>%
            dplyr::select( "MTP", "target.power", 
                           "pt", "w", "power" )
    }
    fin_pts
}



#' @title Return type of pump object (result function)
#' @description Returns whether call was power, mdes, or sample.
#' @return pump_type: power, mdes, or sample, as a string.
#'
#' @rdname pumpresult
#' @export
#' 
pump_type <- function(x)
{
    stopifnot( is.pumpresult(x) || is.pumpgridresult(x) )
    return( attr(x, "type" ) )
}



is_long_table <- function(power_table)
{
    stopifnot( is.pumpresult(power_table) || is.pumpgridresult(power_table) )
    
    lt <- attr( power_table, "long.table" )
    if ( !is.null( lt ) && lt == TRUE ) {
        return( TRUE )
    } else {
        return( FALSE )
    }
}




#' @title Convert power table from wide to long (result function)
#'
#' @description Transform table returned from pump_power 
#' to a long format table or to a wide format table.
#'
#' @param power_table pumpresult object for a power result 
#'  (not mdes or sample). (It can also take a raw dataframe of the wide table 
#'  to convert to long, as an internal helper method.)
#' @param M scalar; set if power_table is a data.frame 
#' without set number of outcomes. Usually ignore this.
#'   
#' @return data.frame of power results in long format.
#'
#' @export
#' 
transpose_power_table <- function(power_table, M = NULL)
{
    
    ptorig <- power_table
    
    pp <- NA
    
    pr <- is.pumpresult(power_table) || is.pumpgridresult(power_table)
    if ( pr ) {
        stopifnot( pump_type( power_table ) == "power" )
        M <- params(power_table)$M
    } else {
        stopifnot( !is.null(M) )
    }
    
    
    if ( !pr || !is_long_table(power_table) ) {
        pnames <- get_power_names(M, long = TRUE)
        
        pp <- power_table %>% 
            as.data.frame() %>%
            tidyr::pivot_longer( cols = tidyselect::any_of( names(pnames) ),
                                 names_to = "power",
                                 values_to = "power_val" ) %>%
            tidyr::pivot_wider( names_from = "MTP",
                                values_from = "power_val" ) %>%
            dplyr::mutate( power = pnames[ .data$power] )
    } else {
        pnames <- get_power_names(M)
        pp <- power_table %>%
            dplyr::mutate( power = pnames[ .data$power ] ) %>%
            tidyr::pivot_longer( 
                cols = tidyselect::any_of(c( "None", 
                                             params(power_table)$MTP )),
                names_to = "MTP",
                values_to = "power_val" ) %>%
            tidyr::pivot_wider( names_from = "power",
                                values_from = "power_val" )
    } 
    
    if ( is.pumpresult( ptorig ) || is.pumpgridresult(ptorig) ) {
        att <- attributes(ptorig)
        att["names"] <- NULL
        att["row.names"] <- NULL
        att["long.table"] <- !att[["long.table"]]
        for (i in seq(1, length(att))) {
            attr( pp, names(att)[[i]] ) <- att[[ i ]]
        }    
    }
    
    return( pp )
}


#' @return is.pumpresult: TRUE if object is a pumpresult object.
#'
#' @export
#'
#' @rdname pumpresult
#' 
is.pumpresult <- function(x)
{
    inherits(x, "pumpresult")
}


#' @return `[`: pull out rows and columns of the dataframe.
#'
#' @rdname pumpresult
#' @export
`[.pumpresult` <- function(x, ...)
{
    as.data.frame(x)[...] 
}




#' @return `[[`: pull out single element of dataframe.
#'
#' @rdname pumpresult
#' @export
`[[.pumpresult` <- function(x, ...)
{
    as.data.frame(x)[[...]] 
}



#' @title Dimension of pumpresult object
#' @return dim: Dimension of pumpresult (as matrix)
#'
#' @rdname pumpresult
#' @export
#' 
dim.pumpresult <- function(x, ... )
{
    return( dim( as.data.frame(x) ) )
}


#' @title Pretty print pump result with parameters
#'
#' @description
#' Calls the print_context method with results and control both set to TRUE.
#'
#' @seealso print_context
#'
#' @param object Object to summarize.
#' @param ... Extra options passed to print.pumpresult
#' @rdname pumpresult
#' 
#' @return summary: No return value; prints results.
#' 
#' @export
#'   
summary.pumpresult <- function(object, ...)
{
    print_context( object, insert_results = TRUE, insert_control = TRUE, ... )
}


calc_binomial_SE <- function(prop, tnum)
{
    pp <- (tnum * prop + 2) / (tnum + 4)
    pp * (1 - pp) / sqrt(tnum + 4)
}


#' @title Pretty print pump result
#'
#' @param ... No extra options passed.
#' @param n Number of lines of search path to print, max.
#' @param header FALSE means skip some header info on the result, just print
#'   the data.frame of actual results.
#' @param search FALSE means don't print the search path for a result for
#'   mdes or sample.
#'   
#' @return print: No return value; prints results.  
#'   
#' @rdname pumpresult
#' 
#' @export
#' 
print.pumpresult <- function(x, n = 10,
                             header = TRUE,
                             search = FALSE,
                             ... )
{
    result_type <- attr( x, "type" )
    
    pars = params(x)
    
    if ( header ) {
        scat( "%s result: %s d_m with %d outcomes\n",
              result_type, d_m(x), pars$M )
        
        if ( result_type == "mdes" || result_type == "sample" ) {
            pow_params <- attr( x, "power.params.list" )
            scat( "  target %s power: %.2f\n", pow_params$power.definition,
                  pow_params$target.power )
        }
    }
    
    tnum <- pars$tnum
    
    if ( is.pumpresult(x) ) {
        
        if ( pump_type(x) == "power" ) {
            SEh <- 0.5 + min( abs( 0.5 - x[,-1] ), na.rm = TRUE )
            SEh <- calc_binomial_SE( SEh, tnum )
            SEl <- 0.5 + max( abs( 0.5 - x[,-1] ), na.rm = TRUE )
            SEl <- calc_binomial_SE( SEl, tnum )
            print( as.data.frame( x ), row.names = FALSE )
            
            if ( pars$M > 1 ) {
                scat("\t%.3f <= SE <= %.3f\n", SEl, SEh )
            }
        } else if ( pump_type(x) == "sample" ) {
            if ( !exact_calc(x) ) {
                SE <- pmax( x[1,4], 1 - x[1,4] )
                SE <- calc_binomial_SE( SE, tnum )
                x$SE <- round( SE, digits = 2 )
            }
            print( as.data.frame( x ), row.names = FALSE )
            
        } else {
            if ( pars$MTP != "None" ) {
                SE <- calc_binomial_SE( x[[3]], tnum )
                x$SE <- round( SE, digits = 2 )
            }
            print( as.data.frame( x ), row.names = FALSE )
        }
    } else {
        nc <- ncol(x)
        pvs <-  x[,nc][[1]]
        x$SE <- calc_binomial_SE( pvs, tnum )
        print( as.data.frame(x) )
    }
    
    if ( search ) {
        print_search( x, n = n )
    } else {
        s_path <- attr( x, "tries" )
        if ( !is.null( s_path ) ) {
            scat( "\t(%d steps in search)\n", nrow(s_path) - 5 )
        }
    }
    
    if ( !is.pumpgridresult(x) && pump_type(x) != "power" && 
         !is.null(attr(x, "flat" )) && attr(x, "flat") ) {
        scat( "Note: Power curve is relatively flat. \n" )
    }
    invisible( x )
}





#' @title Print the search history of a pump result object (result function)
#'
#' @description For pump_mdes and pump_sample, print the 
#' (abbreviated) search history.
#'
#' @inheritParams print.pumpresult
#' 
#' @return No return value; prints results.
#' 
#' @keywords internal
print_search <- function(x, n = 10)
{
    tr <- search_path( x )
    
    if ( !is.null( tr )  ) {
        cat( "\nSearch history\n")
        nr <- nrow( tr )
        if ( nr >= n ) {
            print( utils::head( tr, max(n/2,1) ) )
            scat( "\t...  %s steps total ...\n", nr )
            if ( n >= 2 ) {
                print( utils::tail( tr, n/2 ) )
            }
        } else {
            print( tr )
        }
        invisible( nr )
    } else {
        invisible( 0 )
    }
}



#' @title Print context (design, model, parameter values) of 
#' pumpresult or pumpgridresult (result function)
#'
#' @description
#'
#' Print out the context (design and model, with parameter values) of 
#' given pump result or pump grid result object.  
#' The "***" denotes varying values in the
#' printout.
#' 
#'
#' @param x A pumpresult object or pumpgridresult object.
#' @param insert_results Include actual results in the printout.
#' @param insert_control Include the optimizer control parameter information.
#' @param ... Extra arguments to pass to print.pumpresult.
#' 
#' @return No return value; prints results.
#'
#' @export
#' 
print_context <- function( 
    x, insert_results = FALSE, insert_control = FALSE, ...  ) 
{
    is_grid <- is.pumpgridresult(x)
    
    reduce_vec <- function(vec) {
        if ( is.null(vec) ) {
            return( NULL )
        }
        if ( is.numeric( vec ) ) {
            if ( all( round( vec ) == vec ) ) {
                vec <- as.character( as.integer(vec) )
            } else {
                vec <- as.character( round( vec, digits = 2 ) )
            }
        }
        if ( length( unique( vec ) ) == 1 ) {
            vec[[1]]
        } else {
            paste( vec, collapse = " / " )
        }
    }
    
    params <- params(x)
    MDESv <- params$MDES
    params <- lapply( params, reduce_vec )
    
    d_m <- d_m(x)
    des <- parse_d_m(d_m)
    if ( des$levels < 3 ) {
        params$K <- "none"
    }
    
    result_type <- attr( x, "type" )
    
    if ( result_type == "sample" ) {
        smp_type <- attr(x, "sample.level" )
        if ( smp_type == "nbar" ) {
            params$nbar <- "*"
        } else if ( smp_type == "J" ) {
            params$J <- "*"
        } else if ( smp_type == "K" ) {
            params$K <- "*"
        }
    }
    
    # adjust for grid printout
    if ( is_grid ) {
        print_grid_header( x )
        
        varnames <- attr( x, "var_names" )
        for (v in varnames) {
            params[[v]] <- "***"
        }
        
        if ( "MDES" %in% varnames ) {
            MDESv <- "MDES varying"
        } 
    } else {
        
        scat( "%s result: %s d_m with %s outcomes",
              result_type, d_m(x), params$M )
        if ( !is.null( params$numZero ) ) {
            scat( " (%s zeros)\n", params$numZero )
        } else {
            scat( "\n" )
        }
    }
    
    if ( result_type == "mdes" || result_type == "sample" ) {
        pow_params <- attr( x, "power.params.list" )
        scat( "  target %s power: %.2f\n", pow_params$power.definition,
              pow_params$target.power )
    }
    
    if ( !is.null( MDESv ) ) {
        scat( "  MDES vector: %s\n", paste( MDESv, collapse = ", " ) )
    }
    
    cname <- function(numCov) {
        if ( numCov == "1" ) {
            return( "1 covariate" )
        } else {
            return( sprintf( "%s covariates", numCov ) )
        }
    }
    
    if ( des$levels < 2 ) {
        scat( "  nbar: %s\tTbar: %s\n",
              params$nbar, params$Tbar )
    } else if ( des$levels < 3 ) {
        scat( "  nbar: %s\tJ: %s\tTbar: %s\n",
              params$nbar, params$J, params$Tbar )
    } else {
        scat( "  nbar: %s\tJ: %s\tK: %s\tTbar: %s\n",
              params$nbar, params$J, params$K, params$Tbar )
    }
    scat( "  alpha: %s\t\n", params$alpha)
    scat( "  Level:\n    1: R2: %s (%s)\n",
          params$R2.1, cname( params$numCovar.1 ) )
    
    if ( des$levels >= 2 ) {
        scat( "    2: ")
        if ( des$FE.2 ) {
            scat( "  fixed effects  " )
        } else {
            scat( "R2: %s (%s)", params$R2.2, cname( params$numCovar.2) )
        }
        scat( "\tICC: %s\tomega: %s\n", params$ICC.2, params$omega.2 )
    }
    if ( des$levels >= 3 ) {
        scat( "    3: ")
        if ( des$FE.3 ) {
            scat( "  fixed effects  " )
        } else {
            scat( "R2: %s (%s)", params$R2.3, cname( params$numCovar.3) )
        }
        scat( "\tICC: %s\tomega: %s\n", params$ICC.3, params$omega.3 )
    }
    if ( !is.null( params$rho.matrix ) ) {
        cat( "Rho matrix:\n" )
        print( params$rho.matrix )
    } else {
        scat( "  rho = %s\n", params$rho )
    }
    
    if ( insert_results ) {
        print.pumpresult(x, header = FALSE, search = FALSE, ... )
    }
    
    if ( insert_control ) {
        ex_params <- params[ c("B", "max.steps", 
                               "tnum", "start.tnum", "final.tnum", "tol") ]
        if ( params$MTP != "WY-SS" && params$MTP != "WY-SD" ) {
            ex_params$B <- NULL
        }
        no_inc <- vapply( ex_params, is.null, is.numeric(1) )
        ex_params <- ex_params[ !no_inc ]
        scat( "\t(%s)\n",
              paste( names(ex_params), ex_params, sep = " = ",
                     collapse = "  " ) )
    }
    
    invisible( x )
}





#' @title Cast pumpresult result to data.frame
#'
#' @param row.names NULL or a character vector giving the 
#' row names for the data frame.
#' @param optional logical. If TRUE, setting row names and 
#' converting column names is optional.
#' @param ... additional arguments to be passed to the 
#' as.data.frame.list methods.
#'
#' @return as.data.frame: pumpresult object as a clean 
#' dataframe (no more attributes from pumpresult).
#' @rdname pumpresult
#'
#' @export
#' 
as.data.frame.pumpresult <- function( 
        x, row.names = NULL, optional = FALSE, ... 
) {
    class(x) <- "list"
    as.data.frame( x, row.names = row.names, optional = optional, ... )
    
}


