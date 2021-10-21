



scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#' Update a pump call, tweaking some parameters
#' 
#' @param x Pump result object.
#' @return New call using parameters of old object.
#' 
#' @export
update.pumpresult = function( x, ... ) {
    params = attr(x,"param")
    dts = list(...)
    for ( d in names(dts) ) {
        params[[d]] = dts[[d]]
    }
    params["design"] = design(x)
    result_type = attr( x, "type" )
    if ( result_type == "power" ) {
        do.call(pump_power, params)
    } else if ( result_type == "mdes" ) {
        do.call( pump_mdes, params )
    } else if ( result_type == "sample" ) {
        do.call( pump_sample, params )
    } else {
        stop( sprintf( "Unrecognized type, %s, in update()", result_type ) )
    }
}


make.pumpresult = function( x,
                 type = c( "power", "mdes", "sample" ),
                 design = design,
                 params.list = NULL,
                 tries = NULL, final.pts = NULL,
                 just.result.table = TRUE,
                 ... ) {
    type = match.arg(type)
    class(x) = c( "pumpresult", class(x) )
    attr(x, "type" ) = type
    attr(x, "params.list") = params.list
    attr(x, "design") = design
    ll = list(...)
    for ( l in names(ll) ) {
        attr(x, l) = ll[[ l ]]
    }
    if ( !just.result.table && !is.null( tries ) ) {
        attr( x, "tries" ) = tries
    }
    if ( !is.null( final.pts ) ) {
        attr( x, "final.pts" ) = final.pts
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
    rs$delta = rs$power - rs$target.power
    return( rs )
}


#' Obtain recheck path (to see rate of power change)
#'
#' @param x A pumpresult object
#'
#' @param Dataframe describing power curve.
#' @family pumpresult
#' @export
power_curve = function( x, ... ) {
    stopifnot( is.pumpresult( x ) )
    return( attr( x, "final.pts" ) )
}



#' Dimension of pumpresult
#'
#' @family pumpresult
#' @export
dim.pumpresult = function( x, ... ) {
    c( length( x[[1]] ), length(x) )
}


#' Pretty print pump result
#'
#' @export
#' @param x A pumpresult object.
#' @param ... No extra options passed.
#' @family pumpresult
print.pumpresult = function( x, n = 10, ... ) {
    result_type = attr( x, "type" )

    scat( "%s result: %s design with %d outcomes\n", 
          result_type, design(x), params(x)$M )
    print( as.data.frame( x ) )

    tr = attr( x, "tries" )
    if ( !is.null( tr ) ) {
        cat( "\nSearch history\n")
        nr = nrow( tr )
        if ( nr <= n ) {
            print( tr )
        } else {
            print( head( tr, max(n/2,1) ) )
            scat( "\t...  %s steps total ...\n", nr )
            print( tail( tr, max(n/2),1) )
        }
    }


    invisible( x )
}



#' Print design of given pump result object
#' 
#' @param x A pumpresult object.
#' 
#' @export
print_design = function( x ) {
    
    
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
    params = sapply( params, reduce_vec, simplify=FALSE )
    attach( params, warn.conflicts = FALSE)
    design = design(x)
    result_type = attr( x, "type" )
    des = parse_design(design)
    if ( des$levels < 3 ) {
        params$K = "none"
    }
    
    scat( "%s result: %s design with %s outcomes",
          result_type, design(x), M )
    if ( !is.null( numZero ) ) {
        scat( " (%s zeros)\n", numZero )
    } else {
        scat( "\n" )
    }
    
    if ( !is.null( MDESv ) ) {
        scat( "  MDES vector: %s\n", paste( MDESv, collapse=", " ) )
    }
    
    scat( "  nbar: %s\tJ: %s\tK: %s\tTbar: %s\n", nbar, J, params$K, Tbar )
    scat( "  alpha: %s\t\n", alpha)
    scat( "  Level:\n    1: R2: %s (%s covariate)\n",
          R2.1, numCovar.1 )
    if ( des$levels >= 2 ) {
        scat( "    2: R2: %s (%s covariate)\tICC: %s\tomega: %s\n",
           R2.2, numCovar.2, ICC.2, omega.2 )
    }
    if ( des$levels >= 3 ) {
        scat( "    3: R2: %s (%s covariate)\tICC: %s\tomega: %s\n",
               R2.3, numCovar.3, ICC.3, omega.3 )
    }
    if ( !is.null( rho.matrix ) ) {
        cat( "Rho matrix:\n" )
        print( rho.matrix )
    } else {
        scat( "  rho = %s\n", rho )
    }
    scat( "\t*  B = %s", B )
    if ( exists( "tnum" ) ) {
        scat( "  tnum = %s", tnum)
    } 
    scat( "\n" )
    
    detach( params )
    invisible( params )
}


summary.pumpresult = function( x ) {
    result_type = attr( x, "type" )
    
    scat( "%s result: %s design with %d outcomes\n", 
          result_type, design(x), params(x)$M )
    print( as.data.frame( x ) )
    
    tr = attr( x, "tries" )
    if ( !is.null( tr ) ) {
        cat( "\nSearch history\n")
        nr = nrow( tr )
        if ( nr <= n ) {
            print( tr )
        } else {
            print( head( tr, max(n/2,1) ) )
            scat( "\t...  %s steps total ...\n", nr )
            print( tail( tr, max(n/2),1) )
        }
    }
    
    
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

