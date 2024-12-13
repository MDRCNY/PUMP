#' Parse the power definition
#'
#' @param power.definition i.e. D1indiv, min1, complete
#' @param M number of outcomes
#' @return information about power type
#' @keywords internal
parse_power_definition <- function(power.definition, M) {
    powertype <- list( min = FALSE,
                       complete = FALSE,
                       indiv = FALSE )
    if ( !is.null(power.definition) ) {
        if ( stringr::str_detect( power.definition, "min" ) ) {
            powertype$min <- TRUE
            powertype$min_k <- readr::parse_number( power.definition )
            stopifnot( is.numeric( powertype$min_k ) )
        } else if ( stringr::str_detect( power.definition, "complete" ) ) {
            powertype$min <- TRUE
            powertype$complete <- TRUE
            powertype$min_k <- M
        } else if ( stringr::str_detect( power.definition, "indiv.mean" ) ) {
            powertype$indiv <- TRUE
            powertype$indiv_k <- NULL
        } else if ( stringr::str_detect( power.definition, "indiv" ) ) {
            powertype$indiv <- TRUE
            indiv_k <- readr::parse_number( power.definition )
            if ( !is.na( indiv_k ) ) {
                powertype$indiv_k = indiv_k
            } else {
                stop(glue::glue( "Invalid power definition 
                                 {power.definition}.",
                                 "Try, e.g., 'D1indiv'." ) )
            }
        }
    }
    
    return( powertype )
}



get_power_names <- function(M, long = FALSE) {
    
    if ( M == 1 ) {
        nms <- c( "D1indiv" )
        lnms <- c( "individual outcome 1" )
    } else {
        nms <- c( paste('D', 1:M, "indiv", sep = "" ),
                  'indiv.mean',
                  paste('min', 1:(M - 1), sep = ''),
                  'complete' )
        
        lnms <- c( paste("individual outcome", 1:M),
                   'mean individual',  paste(1:(M - 1),'minimum', sep = '-'),
                   'complete')
    }
    
    if ( long ) {
        names( lnms ) <- nms
        return( lnms )
    } else {
        names( nms ) <- lnms
        return( nms )
    }
}

# # Stolen from development purrr
# silently <- function(.f, otherwise = NULL) {
#     .f <- purrr::as_mapper(.f)
#     function(...) {
#         ret <-
#             purrr:::capture_output(
#                 purrr:::capture_error(.f(...), otherwise, quiet=TRUE)
#             )
#         # reformat output to an un-nested list
#         list(
#             result = ret$result$result,
#             output = ret$output,
#             messages = ret$messages,
#             warnings = ret$warnings,
#             error = ret$result$error
#         )
#     }
# }

# print out results cleanly
scat <- function(str, ...) {
    cat( sprintf( str, ... ) )
}


smessage <- function(str, ...) {
    message( sprintf( str, ... ) )
}


swarning <- function(str, ...) {
    warning( sprintf( str, ... ), call. = FALSE )
}


sstop <- function(str, ...) {
    stop( sprintf( str, ... ), call. = FALSE )
}



