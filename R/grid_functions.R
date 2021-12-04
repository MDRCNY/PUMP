#' Run grid across any of the core pump functions
#'
#' @param args list of scenario arguments
#' @param pum_function pump_mdes, pump_sample, pump_power
#' @param verbose print out detailed diagnostics
#' @param drop_unique_columns TODO
#' @param use_furrr not currently in use, whether to parallelize
#' @param ... Extra arguments passed to the underlying pump_power, pump_sample,
#'   or pump_mdes functions.
run_grid <- function( args, pum_function, verbose = FALSE,
                      drop_unique_columns, ...,
                      use_furrr = FALSE ) {
  grid <- do.call( tidyr::expand_grid, args )
  if ( verbose ) {
    scat( "Processing %d calls\n", nrow(grid) )
  }


  if ( use_furrr ) {
    # TODO: To be fixed later
    grid$res <- furrr::future_pmap( grid, pum_function, ...,
                                   .progress = verbose )
  } else {
    grid$res <- purrr::pmap( grid, pum_function, ...,
                             verbose = verbose )

  }

  if ( drop_unique_columns ) {
    grid <- dplyr::select(
      grid, dplyr::any_of( c("MDES", "numZero" )) |
      tidyselect::vars_select_helpers$where( ~ !is.numeric( .x ) || length( unique( .x ) ) > 1 ) )
  }
  grid$res <- purrr::map( grid$res, as.data.frame )
  grid$MTP <- NULL
  grid <- tidyr::unnest( grid, .data$res )
  if ( "MTP" %in% names(grid) ) {
    grid <- grid %>% dplyr::arrange( .data$MTP ) %>%
    dplyr::relocate( .data$MTP )
  }

  return( grid )
}


#' Setup parallel processing
#'
#' Set up furrr to use all but one core
#'
# @importFrom future plan
# @importFrom future multisession
#' @importFrom parallel detectCores
#'
#' @export
setup_default_parallel_plan <- function() {
  future::plan(future::multisession, workers = parallel::detectCores() - 1 )
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
#' These calls use furrr's future_pmap package to allow for parallel
#' computation.  You can use the `setup_default_parallel_plan()` method or your
#' own.  If you do nothing, it will default to single session.
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
#' @param use_furrr Use parallel processing furrr package to fit grid.
#' @param ... Extra arguments passed to the underlying pump_power, pump_sample,
#'   or pump_mdes functions.
#'
#' @importFrom magrittr %>%
# @importFrom furrr future_pmap
#' @importFrom tidyselect vars_select_helpers
#' @family grid functions
#' @export
pump_power_grid <- function( design, MTP, MDES, M, nbar,
                             J = 1, K = 1, numZero = NULL,
                             Tbar, alpha = 0.05,
                             numCovar.1 = NULL,
                             numCovar.2 = NULL,
                             numCovar.3 = NULL,
                             R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                             ICC.2 = NULL, ICC.3 = NULL,
                             omega.2 = NULL, omega.3 = NULL,
                             rho,
                             verbose = FALSE,
                             use_furrr = FALSE,
                             drop_unique_columns = TRUE,
                             ... ) {

  if ( sum( duplicated( MDES ) ) > 0 ) {
    stop(paste("Cannot pass duplicate MDES values to pump_power_grid.\n",
         "Did you try to give a vector of varying MDES?" ))
  }

  args <- list( design = design, M = M, MDES = MDES,
                J = J, K = K,
                nbar = nbar, numZero = numZero,
                Tbar = Tbar, alpha = alpha,
                numCovar.1 = numCovar.1,
                numCovar.2 = numCovar.2,
                numCovar.3 = numCovar.3,
                R2.1 = R2.1, R2.2 = R2.2,
                ICC.2 = ICC.2, ICC.3 = ICC.3,
                rho = rho,
                omega.2 = omega.2, omega.3 = omega.3 )
  nulls <- purrr::map_lgl( args, is.null )
  args <- args[ !nulls ]

  grid <- run_grid( args, pum_function = pump_power, verbose = verbose,
                    drop_unique_columns = drop_unique_columns,
                    MTP = MTP, ..., use_furrr = use_furrr )

  return( grid )
}



#' Run pump_mdes on combination of factors
#'
#' @inheritParams pump_mdes
#' @inheritParams pump_power_grid
#'
#' @family grid functions
#'
#' @export
pump_mdes_grid <- function( design, MTP, M,
                            target.power, power.definition, tol,
                            nbar, J = 1, K = 1,
                            Tbar, alpha,
                            numCovar.1 = NULL,
                            numCovar.2 = NULL,
                            numCovar.3 = NULL,
                            R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                            ICC.2 = NULL, ICC.3 = NULL,
                            omega.2 = NULL, omega.3 = NULL,
                            rho,
                            verbose = FALSE,
                            use_furrr = FALSE,
                            drop_unique_columns = TRUE,
                            ... ) {


  args <- list( design = design, M = M,
                J = J, K = K, nbar = nbar,
                target.power = target.power,
                power.definition = power.definition,
                MTP = MTP,
                Tbar = Tbar, alpha = alpha,
                numCovar.1 = numCovar.1,
                numCovar.2 = numCovar.2,
                numCovar.3 = numCovar.3,
                R2.1 = R2.1, R2.2 = R2.2,
                ICC.2 = ICC.2, ICC.3 = ICC.3,
                rho = rho,
                omega.2 = omega.2, omega.3 = omega.3 )
  nulls <- purrr::map_lgl( args, is.null )
  args <- args[ !nulls ]

  grid <- run_grid( args, pum_function = pump_mdes,
                    verbose = verbose,
                    drop_unique_columns = drop_unique_columns,
                    tol = tol, ..., use_furrr = use_furrr )

  return( grid )
}


#' Run pump_sample on combination of factors
#'
#' @inheritParams pump_sample
#' @inheritParams pump_power_grid
#'
#' @family grid functions
#'
#' @export
pump_sample_grid <- function( design, MTP, M,
                              target.power, power.definition, tol,
                              MDES = NULL,
                              typesample,
                              nbar = NULL, J = NULL, K = NULL,
                              Tbar, alpha,
                              numCovar.1 = NULL,
                              numCovar.2 = NULL,
                              numCovar.3 = NULL,
                              R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                              ICC.2 = NULL, ICC.3 = NULL,
                              omega.2 = NULL, omega.3 = NULL,
                              rho,
                              verbose = FALSE,
                              use_furrr = FALSE,
                              drop_unique_columns = TRUE,
                              ... ) {


  args <- list( design = design, M = M, J = J, K = K,
                power.definition = power.definition,
                MTP = MTP,
                MDES = MDES, nbar = nbar,
                target.power = target.power,
                Tbar = Tbar, alpha = alpha,
                numCovar.1 = numCovar.1,
                numCovar.2 = numCovar.2,
                numCovar.3 = numCovar.3,
                R2.1 = R2.1, R2.2 = R2.2,
                ICC.2 = ICC.2, ICC.3 = ICC.3,
                rho = rho,
                omega.2 = omega.2, omega.3 = omega.3 )
  nulls <- purrr::map_lgl( args, is.null )
  args <- args[ !nulls ]

  grid <- run_grid( args, pum_function = pump_sample,
                   typesample = typesample,
                   verbose = verbose,
                   drop_unique_columns = drop_unique_columns,
                   tol = tol, ..., use_furrr = use_furrr )

  return( grid )
}
