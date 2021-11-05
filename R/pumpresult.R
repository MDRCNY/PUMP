



scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#' Update a pump call, tweaking some parameters
#'
#' @param x Pump result object.
#' @return New call using parameters of old object.
#'
#' @export
update.pumpresult = function( x, type=NULL, ... ) {
    params = params(x)
    params["design"] = design(x)
    dts = list(...)
    for ( d in names(dts) ) {
        params[[d]] = dts[[d]]
    }

    if ( !is.null( type ) ) {
        result_type = type
        old_type = attr(x, "type" )

        if ( old_type == "sample" ) {
             ss = smp$`Sample size`
             slvl = attr(x, "sample.level" )
             params[[slvl]] = ss
        }

    } else {
        result_type = attr( x, "type" )
    }

    if ( result_type == "power" ) {
        params["target.power"] = NULL
        params["power.definition"] = NULL
        params["tol"] = NULL
        do.call(pump_power, params)
    } else if ( result_type == "mdes" ) {
        params["MDES"] = NULL
        do.call( pump_mdes, params )
    } else if ( result_type == "sample" ) {
        params[params$typesample] = NULL
        do.call( pump_sample, params )
    } else {
        stop( sprintf( "Unrecognized type, %s, in update()", result_type ) )
    }
}



#' Update a pump call to a grid call
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
        do.call( pump_sample_grid, params )
    } else {
        stop( sprintf( "Unrecognized type, %s, in update_grid()", result_type ) )
    }
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
        attr( x, "search.range" ) <- c( min=min( tries$pt, na.rm=TRUE ),
                                        final=tries$pt[ nrow(tries) ],
                                        max=max( tries$pt, na.rm=TRUE ) )
        attr( x, "flat" ) <- flat
    }
    return( x )
}


#' Get parameters for pump result
#'
#' @return List of design parameters used.
#'
#' @family pumpresult
#' @export
params = function( x, ... ) {
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
#' @return design used
#'
#' @family pumpresult
#' @export
design = function( x, ... ) {
    stopifnot( is.pumpresult( x ) )

    pp = attr( x, "design" )
    return( pp )
}


#' Obtain search path of pump_mdes or pump_sample call
#'
#' @param x A pumpresult object
#'
#' @param Dataframe describing search path, if it was saved in the pumpresult object.
#' @family pumpresult
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
#' This is used to see rate of power change.
#'
#' @param x A pumpresult object.
#' @param all Merge in the search path from the original search.
#'
#' @inheritParams estimate_power_curve
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

#' What type of pumpresult
#'
#' @return power, mdes, or sample, as a string.
#'
#' @family pumpresult
#' @export
pump_type = function( x ) {
    return( attr(x, "type" ) )
}

#' Dimension of pumpresult
#'
#' @family pumpresult
#' @export
dim.pumpresult = function( x, ... ) {
    c( length( x[[1]] ), length(x) )
}


#' Pretty print pump result with parameters
#'
#' @export
#' @param x A pumpresult object.
#' @param ... Extra options passed to print.pumpresult
#' @family pumpresult
summary.pumpresult = function( x, ... ) {
    print_design( x, insert_results = TRUE, ... )
}



#' Pretty print pump result
#'
#' @export
#' @param x A pumpresult object.
#' @param ... No extra options passed.
#' @family pumpresult
print.pumpresult = function( x, n = 10, no_header=FALSE, ... ) {
    result_type = attr( x, "type" )

    if ( !no_header ) {
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
    if ( !is.null( tr ) ) {
        cat( "\nSearch history\n")
        nr = nrow( tr )
        if ( nr <= n ) {
            print( tr )
        } else {
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

    scat( "  nbar: %s\tJ: %s\tK: %s\tTbar: %s\n", params$nbar, params$J, params$K, params$Tbar )
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
        print.pumpresult(x, n = n, no_header=TRUE, ... )
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



#' Is object a pumpresult object?
#'
#' @export
#' @aliases pumpresult
#' @param x the object to check.
#' @family pumpresult
is.pumpresult = function( x ) {
    inherits(x, "pumpresult")
}



#' Cast pumpresult result to data.frame
#'
#' @export
#' @aliases pumpresult
#' @param x the pumpresult object to covert
#' @family pumpresult
#'
#' @export
as.data.frame.pumpresult = function( x ) {
    class(x) = "list"
    as.data.frame( x )
}

