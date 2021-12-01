



scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}




make.pumpresult = function( x,
                            type = c( "power", "mdes", "sample" ),
                            design = design,
                            params.list = NULL,
                            tries = NULL,
                            flat = FALSE,
                            ... ) {
    type <- match.arg(type)
    class(x) <- c( "pumpresult", class(x) )
    attr(x, "type" ) <- type
    attr(x, "params.list") <- params.list
    attr(x, "design") <- design
    ll = list(...)
    for ( l in names(ll) ) {
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




#' Update a single pump call to a grid call
#'
#' Given a few lists of parameters, take a pumpresult and call a grid to explore
#' various versions of the initial scenario.
#'
#' @param x Pump result object.
#' @param ... List of parameters to expand into a grid.
#' @return result of calling corresponding grid
#'
#' @export
update_grid = function( x, ... ) {
    params = attr(x,"param")
    params["design"] = design(x)
    for ( p in names(params) ) {
        params[[p]] = unique( params[[p]] )
    }
    pparam = attr( x, "power.params.list" )
    params = c( params, pparam )

    dts = list(...)
    for ( d in names(dts) ) {
        params[[d]] = dts[[d]]
    }
    result_type = attr( x, "type" )
    if ( result_type == "power" ) {
        params["MDES"] = unique(params[["MDES"]])
        do.call(pump_power_grid, params)
    } else if ( result_type == "mdes" ) {
        do.call( pump_mdes_grid, params )
    } else if ( result_type == "sample" ) {
        params["MDES"] = unique(params[["MDES"]])
        if ( is.null( params[["typesample"]] ) ) {
            params["typesample"] = attr( x, "sample.level" )
        }
        params[params$typesample] = NULL
        do.call( pump_sample_grid, params )
    } else {
        stop( sprintf( "Unrecognized type, %s, in update_grid()", result_type ) )
    }
}






#' Update a pump call, tweaking some parameters
#'
#' One of the optional parameters can be a `type = something` argument, where
#' the "something" is either "power", "sample", or "mdes", if the call should be
#' shifted to a different pump call (pump_power, pump_sample, or pump_mdes,
#' respectively).
#'
#' @param object Pump result object.
#' @param ... Parameters as specified in `pump_power`, `pump_mdes`, and
#'   `pump_sample` that should be overwritten.
#'
#' @return Results of a new call using parameters of old object with newly
#'   specified parameters replaced.
#'
#' @export
update.pumpresult = function( object, ... ) {
    params = params(object)
    params["design"] = design(object)
    result_type = attr(object, "type" )
    params["type"] = result_type

    dts = list(...)
    for ( d in names(dts) ) {
        params[[d]] = dts[[d]]
    }

    # Are we changing what kind of calculation we want to perform?  If so,
    # adjust some parameters as needed.
    if ( params$type != result_type ) {

        if ( result_type == "sample" ) {
            ss = object$`Sample.size`
            slvl = attr(object, "sample.level" )
            params[[slvl]] = ss
        }
        result_type = params$type
    }
    params$type = NULL

    if ( result_type == "power" ) {
        params["target.power"] = NULL
        params["power.definition"] = NULL
        params["tol"] = NULL
        do.call(pump_power, params)
    } else if ( result_type == "mdes" ) {
        params["MDES"] = NULL
        params["numZero"] = NULL
        do.call( pump_mdes, params )
    } else if ( result_type == "sample" ) {
        if ( is.null( params[["typesample"]] ) ) {
            params["typesample"] = attr( object, "sample.level" )
        }
        params[params$typesample] = NULL
        params["numZero"] = NULL
        do.call( pump_sample, params )
    } else {
        stop( sprintf( "Unrecognized type, %s, in update()", result_type ) )
    }
}






#' @title pumpresult object for results of power calculations
#' @name pumpresult
#'
#' @description
#' The pumpresult object is an S3 class that holds the results from
#' `pump_power()`, `pump_sample()`, and `pump_mdes()`.
#'
#' It has several methods that pull different information from this object, and
#' some printing methods for getting nicely formatted results.
#'
#' Pump result objects are also data.frames, so they can be easily manipulated
#' and combined.  The return values from the `grid` functions will just return
#' data frames in general.
#'
#' @seealso updateå
#' @seealso update_grid
#'
#' @param x a pumpresult object (except for is.pumpresult, where it is a generic
#'   object to check).
#' @rdname pumpresult
NULL




#' Get parameters for pump result
#'
#' @return params: List of design parameters used.
#'
#' @rdname pumpresult
#' @export
params <- function( x, ... ) {
    stopifnot( is.pumpresult( x ) )

    pp = attr( x, "params.list" )
    pp_pow = attr(x, "power.params.list" )
    if ( !is.null( pp_pow ) ) {
        pp = c( pp, pp_pow )
    }
    return( pp )
}




#' Get design for pump result
#'
#' @return design: design used (as string)
#'
#' @rdname pumpresult
#' @export
design <- function( x, ... ) {
    stopifnot( is.pumpresult( x ) )

    pp <- attr( x, "design" )
    return( pp )
}


#' Obtain search path of pump_mdes or pump_sample call
#'
#' @return search_path: Dataframe describing search path, if it was saved in the pumpresult object.
#' @rdname pumpresult
#'
#' @export
search_path = function( x, ... ) {
    stopifnot( is.pumpresult( x ) )
    rs <- attr( x, "tries" )
    if ( !is.null( rs ) ) {
        rs$delta = rs$power - rs$target.power
    }
    return( rs )
}


#' Obtain power curve over a range of parameters
#'
#' This is used to see rate of power change as a function of sample size or
#' MDES.
#'
#' @param x A pumpresult object.
#'
#' @param all If TRUE, merge in the search path from the original search.
#' @param low Low range for the plot x-axis.
#' @param high High range for the plot.
#' @param grid.size Number of points to calculate power.
#' @param tnum Number of iterations to calculate power at each grid point.
#'
#' @importFrom rlang .data
#'
#' @export
power_curve <- function( x, all = FALSE,
                         low = NULL, high = NULL, grid.size = 5, tnum = 2000 ) {
    stopifnot( is.pumpresult( x ) )
    fin_pts <- attr( x, "final.pts" )
    if ( is.null( fin_pts ) ) {
        fin_pts = estimate_power_curve( x )
    }
    srch = search_path( x )
    if ( all && !is.null( srch ) ) {
        srch = dplyr::filter( srch, .data$pt < max( fin_pts$pt * 1.1 ) )
        fin_pts = dplyr::bind_rows( fin_pts, srch ) %>%
            dplyr::arrange( .data$pt ) %>%
            dplyr::select( .data$MTP, .data$target.power, .data$pt, .data$w, .data$power )
    }
    fin_pts
}

#' @return pump_type: power, mdes, or sample, as a string.
#'
#' @rdname pumpresult
#' @export
pump_type <- function( x ) {
    return( attr(x, "type" ) )
}



#' @return is.pumpresult: TRUE if object is a pumpresult object.
#'
#' @export
#'
#' @rdname pumpresult
is.pumpresult = function( x ) {
    inherits(x, "pumpresult")
}




#' @return dim: Dimension of pumpresult (as matrix)
#'
#' @rdname pumpresult
#' @export
dim.pumpresult <- function( x, ... ) {
    c( length( x[[1]] ), length(x) )
}


#' Pretty print pump result with parameters
#'
#' @export
#' @param object Object to summarize.
#' @param ... Extra options passed to print.pumpresult
#' @rdname pumpresult
summary.pumpresult = function( object, ... ) {
    print_design( object, insert_results = TRUE, ... )
}



#' Pretty print pump result
#'
#' @export
#' @param ... No extra options passed.
#' @param n Number of lines of search path to print, max.
#' @param header FALSE means skip some header info on the result, just print
#'   the data.frame of actual results.
#' @param search FALSE means don't print the search path for a result for
#'   mdes or sample.
#' @rdname pumpresult
print.pumpresult = function( x, n = 10,
                             header=TRUE,
                             search = header,
                             ... ) {
    result_type = attr( x, "type" )

    if ( header ) {
        scat( "%s result: %s design with %d outcomes\n",
              result_type, design(x), params(x)$M )

        if ( result_type == "mdes" || result_type == "sample" ) {
            pow_params = attr( x, "power.params.list" )
            scat( "  target %s power: %.2f\n", pow_params$power.definition,
                  pow_params$target.power )
        }
    }

    print( as.data.frame( x ), row.names=FALSE )


    tr = attr( x, "tries" )
    if ( search && !is.null( tr ) ) {
        cat( "\nSearch history\n")
        nr = nrow( tr )
        if ( nr <= n ) {
            print( tr )
            print( utils::head( tr, max(n/2,1) ) )
            scat( "\t...  %s steps total ...\n", nr )
            print( utils::tail( tr, max(n/2),1) )
        }
    }


    invisible( x )
}



#' Print design of given pump result object
#'
#' @param x A pumpresult object.
#' @param insert_results Include actual results in the printout.
#' @param ... Extra arguments to pass to print.pumpresult.
#'
#' @export
print_design <- function( x, insert_results = FALSE, ...  ) {


    reduce_vec = function( vec ) {
        if ( is.null(vec) ) {
            return( NULL )
        }
        if ( is.numeric( vec ) ) {
        if ( all( round( vec ) == vec ) ) {
            vec = as.character( as.integer(vec) )
        } else {
            vec = as.character( round( vec, digits=2 ) )
        }
        }
        if ( length( unique( vec ) ) == 1 ) {
            vec[[1]]
        } else {
            paste( vec, collapse=" / " )
        }
    }

    params = params(x)
    MDESv = params$MDES
    params = sapply( params, reduce_vec, simplify = FALSE )

    design = design(x)
    des = parse_design(design)
    if ( des$levels < 3 ) {
        params$K = "none"
    }

    result_type = attr( x, "type" )

    if ( result_type == "sample" ) {
        smp_type = attr(x, "sample.level" )
        if ( smp_type == "nbar" ) {
            params$nbar = "*"
        } else if ( smp_type == "J" ) {
            params$J = "*"
        } else if ( smp_type == "K" ) {
            params$K = "*"
        }
    }

    scat( "%s result: %s design with %s outcomes",
          result_type, design(x), params$M )
    if ( !is.null( params$numZero ) ) {
        scat( " (%s zeros)\n", params$numZero )
    } else {
        scat( "\n" )
    }

    if ( result_type == "mdes" || result_type == "sample" ) {
        pow_params = attr( x, "power.params.list" )
        scat( "  target %s power: %.2f\n", pow_params$power.definition,
              pow_params$target.power )
    }

    if ( !is.null( MDESv ) ) {
        scat( "  MDES vector: %s\n", paste( MDESv, collapse=", " ) )
    }

    scat( "  nbar: %s\tJ: %s\tK: %s\tTbar: %s\n",
          params$nbar, params$J, params$K, params$Tbar )
    scat( "  alpha: %s\t\n", params$alpha)
    scat( "  Level:\n    1: R2: %s (%s covariate)\n",
          params$R2.1, params$numCovar.1 )
    if ( des$levels >= 2 ) {
        scat( "    2: ")
        if ( des$FE.2 ) {
            scat( "  fixed effects  " )
        } else {
            scat( "R2: %s (%s covariate)", params$R2.2, params$numCovar.2)
        }
        scat( "\tICC: %s\tomega: %s\n", params$ICC.2, params$omega.2 )
    }
    if ( des$levels >= 3 ) {
        scat( "    3: ")
        if ( des$FE.3 ) {
            scat( "  fixed effects  " )
        } else {
            scat( "R2: %s (%s covariate)", params$R2.3, params$numCovar.3)
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
        print.pumpresult(x, header=FALSE, search=FALSE, ... )
    }

    scat( "\t  (B = %s", params$B )
    if ( exists( "tnum" ) ) {
        scat( "  tnum = %s", params$tnum)
    }
    if ( exists( "tol" ) ) {
        scat( "  tol = %s", params$tol )
    }
    scat( ")\n" )

    invisible( x )
}





#' Cast pumpresult result to data.frame
#'
#' @param row.names NULL or a character vector giving the row names for the data frame.
#' @param optional logical. If TRUE, setting row names and converting column names is optional.
#' @param ... additional arguments to be passed to the as.data.frame.list methods.
#'
#' @return as.data.frame: pumpresult object as a clean dataframe (no more attributes from
#'   pumpresult).
#' @rdname pumpresult
#'
#' @export
as.data.frame.pumpresult = function( x, row.names = NULL, optional = FALSE, ... ) {
    class(x) = "list"
    as.data.frame( x, row.names = row.names, optional = optional, ... )

}


