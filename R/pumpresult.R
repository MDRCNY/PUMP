



scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}




make.pumpresult = function( x,
                 type = c( "power", "mdes", "sample" ),
                 params.list = NULL,
                 tries = NULL, final.pts = NULL,
                 just.result.table = TRUE,
                 ... ) {
    type = match.arg(type)
    class(x) = c( "pumpresult", class(x) )
    attr(x, "type" ) = type
    attr(x, "params.list") = params.list
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


#' Obtain search path of pump_mdes or pump_sample call
#'
#' @param x A pumpresult object
#'
#' @param Dataframe describing search path, if it was saved in the pumpresult object.
#' @family pumpresult
#' @export
search_path = function( x, ... ) {
    stopifnot( is.pumpresult( x ) )
    return( attr( x, "tries" ) )
}

#' Obtain search path of pump_mdes or pump_sample call
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
    args = attr( x, "args" )

    result_type = args$type

    print( as.data.frame( x ) )

    tr = attr( x, "tries" )
    if ( !is.null( tr ) ) {
        cat( "\nSearch history\n")
        nr = nrow( tr )
        if ( nr < n ) {
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

