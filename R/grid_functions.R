#' Run grid across any of the core pump functions
#'
#' @param args list of scenario arguments
#' @param pum_function pump_mdes, pump_sample, pump_power
#' @param verbose print out detailed diagnostics
#' @param drop.unique.columns logical; drop all parameter columns 
#' that did not vary across the grid.
#' @param use.furrr not currently implemented; whether to use furr package
#' for parallelization
#' @param ... Extra arguments passed to the underlying pump_power, pump_sample,
#'   or pump_mdes functions.
#'   
#' @keywords internal
run_grid <- function(args, pum_function, verbose = FALSE,
                     drop.unique.columns, ...,
                     use.furrr = FALSE) {
  
  # check for duplicate values
  lens <- purrr::map_dbl( args, length )
  for (nm in names(lens)[lens > 1]) {
    if ( length( args[[nm]] ) != length( unique( args[[nm]] ) ) ) {
      stop( sprintf(
          "Cannot pass repeats of same value of %s in a grid call.", nm
      ) )
    }
  }
  
  grid <- do.call( tidyr::expand_grid, args )
  if ( verbose ) {
    scat( "Processing %d calls\n", nrow(grid) )
  }

  if ( use.furrr ) {
    # TODO: To be fixed later
    grid$res <- furrr::future_pmap( grid, pum_function, ...,
                                   .progress = verbose )
  } else {
    grid$res <- purrr::pmap( grid, pum_function, ...,
                             verbose = verbose )

  }

  cnts <- dplyr::summarise( grid,
    dplyr::across( dplyr::everything(),  ~ length( unique( .x ) ) ) )
  var_names <- names(cnts)[ as.numeric( cnts ) > 1 ]
  var_names <- setdiff( var_names, "res" )
  
  if ( verbose ) {
    message( sprintf( "Grid search across multiple %s\n",
                      paste0( var_names, collapse = ", " ) ) )
  }
  
  if ( drop.unique.columns ) {
    grid <- dplyr::select( grid,
      dplyr::any_of(
        unique( c("d_m", "MDES", "numZero", "propZero", var_names, "res" ) ) ) 
      )
  }
  
  params <- params( grid$res[[1]] )
  for (v in var_names) {
    params[v] <- "***"
  }
  
  grid$res <- purrr::map( grid$res, as.data.frame )
  grid$MTP <- NULL
  
  corenames <- colnames( grid$res[[1]] )
  
  grid <- tidyr::unnest( grid, "res" )
  if ( "MTP" %in% names(grid) ) {
    grid <- grid %>% dplyr::arrange( "MTP" ) %>%
    dplyr::relocate( "MTP" )
  }

  attr( grid, "var_names" ) <- var_names 
  attr( grid, "params.list" ) <- params
  
  return( grid )
}


#' Setup parallel processing
#'
#' Set up furrr to use all but one core
#'
#' @importFrom future plan multisession
#' @importFrom parallel detectCores
#' @keywords internal
setup_default_parallel_plan <- function() {
  future::plan(future::multisession, workers = parallel::detectCores() - 1 )
}


#' @title Run pump_power on varying values of parameters (grid
#'   function)
#'
#' @description This extension of `pump_power()` will take lists of
#'   parameter values and run `pump_power()` on all combinations of
#'   these values.
#'
#'   It can only assume the same MDES value for all outcomes due to
#'   this.  (I.e., a vector of MDES values will be interpreted as a
#'   sequence of calls to pump_power, one for each MDES value given).
#'
#'   Each parameter in the parameter list can be a list, not scalar.
#'   It will cross all combinations of the list.
#'
#' @inheritParams pump_power
#'
#' @param MDES vector of numeric; This is *not* a list of MDES for
#'   each outcome, but rather a list of MDES to explore. Each value
#'   will be assumed held constant across all M outcomes.
#' @param verbose logical; TRUE means print out some text as calls
#'   processed.  FALSE do not.
#' @param propZero Proportion of outcomes that have 0 impact (this
#'   will be used to override numZero, only one can be defined)
#' @param drop.unique.columns logical; drop all parameter columns that
#'   did not vary across the grid.
#' @param ... extra arguments passed to the underlying pump_power,
#'   pump_sample, or pump_mdes functions.
#'
#' @return a pumpgridresult object containing power results.
#'
#' @importFrom magrittr %>%
#' @family grid functions
#' @export
#'
#' @examples
#' g <- pump_power_grid( d_m = "d3.2_m3ff2rc", MTP = c( "HO", "BF" ),
#'  MDES = 0.10, J = seq(5, 10, 1), M = 5, K = 7, nbar = 58,
#'  Tbar = 0.50, alpha = 0.15, numCovar.1 = 1,
#'  numCovar.2 = 1, R2.1 = 0.1, R2.2 = 0.7,
#'  ICC.2 = 0.25, ICC.3 = 0.25, rho = 0.4, tnum = 1000)
pump_power_grid <- function(d_m, MTP = NULL, MDES, M = 1, nbar,
                            J = 1, K = 1,
                            propZero = NULL, numZero = NULL,
                            Tbar, alpha = 0.05,
                            numCovar.1 = NULL,
                            numCovar.2 = NULL,
                            numCovar.3 = NULL,
                            R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                            ICC.2 = NULL, ICC.3 = NULL,
                            omega.2 = NULL, omega.3 = NULL,
                            rho = NULL,
                            long.table = FALSE,
                            verbose = FALSE,
                            drop.unique.columns = TRUE,
                            ... ) 
{

  if ( sum( duplicated( MDES ) ) > 0 ) {
    stop(paste("Cannot pass duplicate MDES values to pump_power_grid.\n",
         "Did you try to give a vector of varying MDES?" ))
  }
    
  args <- list( d_m = d_m, M = M, MDES = MDES,
                J = J, K = K,
                nbar = nbar, numZero = numZero,
                propZero = propZero,
                Tbar = Tbar, alpha = alpha,
                numCovar.1 = numCovar.1,
                numCovar.2 = numCovar.2,
                numCovar.3 = numCovar.3,
                R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                ICC.2 = ICC.2, ICC.3 = ICC.3,
                rho = rho,
                omega.2 = omega.2, omega.3 = omega.3 )
  nulls <- purrr::map_lgl( args, is.null )
  args <- args[ !nulls ]

  grid <- run_grid( args, pum_function = pump_power,
                    long.table = long.table,
                    verbose = verbose,
                    drop.unique.columns = drop.unique.columns,
                    MTP = MTP, ..., use.furrr = FALSE )

  args$MTP <- MTP
  args <- c( args, list(...) )
  grid <- make.pumpgridresult(
    grid, "power",
    d_m = d_m,
    long.table = long.table )

  return( grid )
}



#' @title Run pump_mdes on varying values of parameters (grid function)
#' 
#' @description See pump_power_grid() for more details.
#'
#' @inheritParams pump_mdes
#' @inheritParams pump_power_grid
#' 
#' @return a pumpgridresult object containing MDES results.
#'
#' @family grid functions
#'
#' @export
#' 
#' @examples
#' g <- pump_mdes_grid(d_m = "d3.2_m3ff2rc", MTP = "HO",
#'   target.power = c( 0.50, 0.80 ), power.definition = "D1indiv",
#'   tol = 0.05, M = 5, J = c( 3, 9 ), K = 7, nbar = 58,
#'   Tbar = 0.50, alpha = 0.15, numCovar.1 = 1, numCovar.2 = 1,
#'   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.9,
#'   rho = 0.4, tnum = 200)
pump_mdes_grid <- function(d_m, MTP = NULL, M = 1,
                           target.power = 0.80, power.definition = NULL, tol = 0.01,
                           propZero = NULL, numZero = NULL,
                           nbar, J = 1, K = 1,
                           Tbar = 0.5, alpha = 0.05,
                           numCovar.1 = NULL,
                           numCovar.2 = NULL,
                           numCovar.3 = NULL,
                           R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                           ICC.2 = NULL, ICC.3 = NULL,
                           omega.2 = NULL, omega.3 = NULL,
                           rho = NULL,
                           verbose = FALSE,
                           drop.unique.columns = TRUE,
                           ...) {


  args <- list( d_m = d_m, M = M,
                J = J, K = K, nbar = nbar,
                target.power = target.power,
                power.definition = power.definition,
                MTP = MTP,
                numZero = numZero,
                Tbar = Tbar, alpha = alpha,
                numCovar.1 = numCovar.1,
                numCovar.2 = numCovar.2,
                numCovar.3 = numCovar.3,
                R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                ICC.2 = ICC.2, ICC.3 = ICC.3,
                rho = rho,
                omega.2 = omega.2, omega.3 = omega.3 )
  nulls <- purrr::map_lgl( args, is.null )
  args <- args[ !nulls ]

  grid <- run_grid( args, pum_function = pump_mdes,
                    verbose = verbose,
                    drop.unique.columns = drop.unique.columns,
                    tol = tol, ..., use.furrr = FALSE )
  
  args <- c( args, list(...) )
  
  grid <- make.pumpgridresult(
      grid, "mdes",
      d_m = d_m )

  return( grid )
}


#' @title Run pump_sample on varying values of parameters (grid function)
#' 
#' @description See pump_power_grid() for further details.
#'
#' @inheritParams pump_sample
#' @inheritParams pump_power_grid
#' 
#' @return a pumpgridresult object containing sample results.
#'
#' @family grid functions
#'
#' @export
#'
#' @examples
#' g <- pump_sample_grid(d_m = "d3.2_m3ff2rc", typesample = "J",
#'   MTP = "HO", MDES = 0.10, target.power = c( 0.50, 0.80 ),
#'   power.definition = "min1", tol = 0.03,
#'   M = 5, K = 7, nbar = 58, Tbar = 0.50,
#'   alpha = 0.15, numCovar.1 = 1, numCovar.2 = 1,
#'   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.25, ICC.3 = 0.25,
#'   rho = 0.4, tnum = 400)
#'  
pump_sample_grid <- function(d_m, MTP = NULL, M = 1,
                             target.power, power.definition, tol = 0.01,
                             MDES = NULL,
                             propZero = NULL, numZero = NULL,
                             typesample,
                             nbar = NULL, J = NULL, K = NULL,
                             Tbar, alpha,
                             numCovar.1 = NULL,
                             numCovar.2 = NULL,
                             numCovar.3 = NULL,
                             R2.1 = NULL, R2.2 = NULL, R2.3 = NULL,
                             ICC.2 = NULL, ICC.3 = NULL,
                             omega.2 = NULL, omega.3 = NULL,
                             rho = NULL,
                             verbose = FALSE,
                             drop.unique.columns = TRUE,
                             ...) 
{
  args <- list( d_m = d_m, M = M, J = J, K = K,
                power.definition = power.definition,
                MTP = MTP,
                MDES = MDES, nbar = nbar,
                target.power = target.power,
                numZero = numZero,
                Tbar = Tbar, alpha = alpha,
                numCovar.1 = numCovar.1,
                numCovar.2 = numCovar.2,
                numCovar.3 = numCovar.3,
                R2.1 = R2.1, R2.2 = R2.2, R2.3 = R2.3,
                ICC.2 = ICC.2, ICC.3 = ICC.3,
                rho = rho,
                omega.2 = omega.2, omega.3 = omega.3 )
  nulls <- purrr::map_lgl( args, is.null )
  args <- args[ !nulls ]

  grid <- run_grid( args, pum_function = pump_sample,
                   typesample = typesample,
                   verbose = verbose,
                   drop.unique.columns = drop.unique.columns,
                   tol = tol, ..., use.furrr = FALSE )

  args <- c( args, list(...) )
  
  grid <- make.pumpgridresult(
      grid, "sample",
      d_m = d_m,
      sample.level = typesample )
  
  return( grid )
}
