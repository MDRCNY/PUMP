
##
## This file has the "grid" functions that call pump_power, pump_mdes, and pump_sample on multiple combinations
##

scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}


run_grid = function( args, method, verbose, drop_unique_columns, ... ) {
  grid = do.call( expand.grid, args )
  if ( verbose ) {
    scat( "Processing %d calls\n", nrow(grid) )
  }

  grid$res = purrr::pmap( grid, method, ... )

  if ( drop_unique_columns ) {
    grid = dplyr::select( grid, where( ~ !is.numeric( .x ) || length( unique( .x ) ) > 1 ) )
  }

  grid$res = purrr::map( grid$res, as.data.frame )
  grid$res = purrr::map( grid$res, tibble::rownames_to_column, var ="adjustment" )
  grid = tidyr::unnest( grid, res ) %>% dplyr::arrange( adjustment ) %>%
    dplyr::relocate( adjustment )

  grid
}



#' Run pump_power on combination of factors
#'
#' This extenstion of `pump_power()` will take lists of parameter values and run
#' `pump_power()` on all combinations of these values.
#'
#' It can only assume the same MDES value for all outcomes due to this.  (I.e.,
#' a vector of MDES values will be interpreted as a sequence of calls to
#' pump_power, one for each MDES value given).
#'
#' Each parameter in the parameter list can be a list, not scalar.  It will
#' cross all combinations of the list.
#'
#' @inheritParams pump_power
#'
#' @param MDES This is *not* a list of MDES for each outcome, but rather a list
#'   of MDES to explore. Each value will be assumed held constant across all M
#'   outcomes.
#' @param verbose TRUE means print out some text as calls processed.  FALSE do
#'   not.
#' @param drop_unique_columns Drop all parameter colunms that did not vary
#'   across the grid.
#'
#' @export
pump_power_grid <- function( design, MTP, MDES, M, nbar, J = 1, K = 1, numZero = NULL, Tbar, alpha,
                             numCovar.1 = NULL, numCovar.2 = NULL, numCovar.3 = NULL,
                             R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                             ICC.2 = NULL, ICC.3 = NULL,
                             omega.2 = NULL, omega.3 = NULL,
                             rho,
                             verbose = FALSE,
                             drop_unique_columns = TRUE,
                             ... ) {


  args = list( M = M, MDES = MDES, J = J, K = K, nbar = nbar, numZero = numZero,
               Tbar = Tbar, alpha = alpha,
               numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
               R2.1 = R2.1, R2.2 = R2.2, ICC.2 = ICC.2, ICC.3 = ICC.3,
               rho = rho, omega.2 = omega.2, omega.3 = omega.3 )
  nulls = purrr::map_lgl( args, is.null )
  args = args[ !nulls ]

  grid = run_grid( args, method=pump_power, verbose=verbose, drop_unique_columns = drop_unique_columns,
                   design = design, MTP = MTP, ... )

  grid
}







#' Run pump_mdes on combination of factors
#'
#' @inheritParams pump_mdes
#'
#' @family pump_power_grid
#'
#' @export
pump_mdes_grid <- function( design, MTP, M,
                            target.power, power.definition, tol,
                            nbar, J = 1, K = 1, numZero = NULL, Tbar, alpha,
                             numCovar.1 = NULL, numCovar.2 = NULL, numCovar.3 = NULL,
                             R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                             ICC.2 = NULL, ICC.3 = NULL,
                             omega.2 = NULL, omega.3 = NULL,
                             rho,
                             verbose = FALSE,
                             drop_unique_columns = TRUE,
                             ... ) {


  args = list( M = M, J = J, K = K, nbar = nbar, target.power = target.power,
               Tbar = Tbar, alpha = alpha, numZero = numZero,
               numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
               R2.1 = R2.1, R2.2 = R2.2, ICC.2 = ICC.2, ICC.3 = ICC.3,
               rho = rho, omega.2 = omega.2, omega.3 = omega.3 )
  nulls = purrr::map_lgl( args, is.null )
  args = args[ !nulls ]

  grid = run_grid( args, method = pump_mdes, power.definition = power.definition,
                   verbose=verbose, drop_unique_columns = drop_unique_columns,
                   tol = tol, design=design, MTP = MTP, just.result.table = TRUE )

  grid
}








#' Run pump_sample on combination of factors
#'
#' @inheritParams pump_mdes
#'
#' @family pump_power_grid
#'
#' @export
pump_sample_grid <- function( design, MTP, M,
                            target.power, power.definition, tol,
                            MDES = NULL,
                            typesample,
                            nbar = NULL, J = NULL, K = NULL, numZero = NULL, Tbar, alpha,
                            numCovar.1 = NULL, numCovar.2 = NULL, numCovar.3 = NULL,
                            R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                            ICC.2 = NULL, ICC.3 = NULL,
                            omega.2 = NULL, omega.3 = NULL,
                            rho,
                            verbose = FALSE,
                            drop_unique_columns = TRUE,
                            ... ) {


  args = list( M = M, J = J, K = K, MDES = MDES, nbar = nbar, target.power = target.power,
               Tbar = Tbar, alpha = alpha, numZero = numZero,
               numCovar.1 = numCovar.1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3,
               R2.1 = R2.1, R2.2 = R2.2, ICC.2 = ICC.2, ICC.3 = ICC.3,
               rho = rho, omega.2 = omega.2, omega.3 = omega.3 )
  nulls = purrr::map_lgl( args, is.null )
  args = args[ !nulls ]

  grid = run_grid( args, method = pump_sample, power.definition = power.definition,
                   typesample = typesample,
                   verbose=verbose, drop_unique_columns = drop_unique_columns,
                   tol = tol, design=design, MTP = MTP, just.result.table = TRUE )

  grid
}




