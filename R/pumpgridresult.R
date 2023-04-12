


make.pumpgridresult <- function(x,
                                type = c( "power", "mdes", "sample" ),
                                d_m = d_m,
                                params.list = NULL,
                                ...) {
    type <- match.arg(type)
    class(x) <- c( "pumpgridresult", class(x) )
    attr(x, "type" ) <- type
    if ( !is.null( params.list ) ) {
        attr(x, "params.list") <- params.list
    } else {
        stopifnot( !is.null( attr(x,"params.list" ) ) )
    }
    attr(x, "d_m") <- d_m

    ll <- list(...)
    for (l in names(ll)) {
        attr(x, l) <- ll[[ l ]]
    }
    
    return( x )
}




#' @title Result object for results of grid power calculations
#' 
#' @name pumpgridresult
#'
#' @description
#' The pumpgridresult object is an S3 class that holds the results from
#' `pump_power_grid()`, `pump_sample_grid()`, and `pump_mdes_grid()`.
#'
#' It has several methods that pull different information from this object, and
#' some printing methods for getting nicely formatted results.
#'
#'
#' @param x a pumpgridresult object 
#' (except for is.pumpgridresult, where it is a generic object to check).
#' @rdname pumpgridresult
NULL






#' @return is.pumpgridresult: TRUE if object is a pumpgridresult object.
#'
#' @rdname pumpgridresult
#' 
#' @export
is.pumpgridresult <- function(x) {
    inherits(x, "pumpgridresult")
}



print_grid_header <- function(x) {
    result_type <- attr( x, "type" )
    
    d_m <- d_m(x)
    if ( length( d_m ) > 1 ) {
        d_m <- paste0( "multi-design ", paste( d_m, collapse = "/" ) )
    }
    scat( "%s grid result: %s d_m with %s outcomes\n",
          result_type, d_m(x), params(x)$M )
    
    scat( "Varying across %s\n",
          paste0( attr( x, "var_names" ), collapse = ", " ) )
    
    if ( result_type == "mdes" || result_type == "sample" ) {
        pow_params <- attr( x, "power.params.list" )
        scat( "  target %s power: %.2f\n", pow_params$power.definition,
              pow_params$target.power )
    }
}


#' @title Pretty print pump grid result (result function)
#'
#' @param ... extra options passed.
#' @param header logical; FALSE means skip some 
#' header info on the result, just print
#' the data.frame of actual results.
#' @rdname pumpgridresult
#' 
#' @return print: No return value; prints results.
#' 
#' @export
print.pumpgridresult <- function(x,header = TRUE, ...) 
{
    if ( header ) {
       print_grid_header( x )
    }
    
    print( as.data.frame( x ), row.names = FALSE )
    
    invisible( x )
}

#' @title Pretty print pump grid result with parameters (result function)
#'
#' @param object object to summarize.
#' @param ... extra options passed to print.pumpgridresult
#' @rdname pumpgridresult
#' 
#' @return summary: No return value; prints results.
#' 
#' @export
summary.pumpgridresult <- function(object, ...)
{
    print_grid_header( object )
    
    print_context( object, 
                  insert_results = TRUE, insert_control = TRUE, ... )
    
    invisible( object )
}