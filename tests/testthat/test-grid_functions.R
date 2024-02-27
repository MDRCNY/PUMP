# library( PUMP )
# library( testthat )

test_that("pump_power_grid works", {

  
  pp <- pump_power_grid(    d_m = "d3.2_m3ff2rc",
                            MTP = "HO",
                            MDES = c( 0.05, 0.2 ),
                            M = 5,
                            J = c( 3, 5, 9), # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200,
                            verbose = FALSE
  )
  pp
  expect_equal( nrow(pp), 3 * 2 * 2 )

})

test_that("More extensive testing of grid()", {
 
  skip_on_cran()
  pp <- pump_power_grid(    d_m = "d3.2_m3ff2rc",
                            MTP = "HO",
                            MDES = c( 0.05, 0.2 ),
                            M = 5,
                            J = 4, # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200,
                            verbose = FALSE
  )
  expect_equal( nrow(pp), 2 * 2 )








  pp <- pump_power_grid(    d_m = "d3.2_m3ff2rc",
                            MTP = "HO",
                            MDES = c( 0.05, 0.2 ),
                            numZero = c(1,2,3),
                            M = 5,
                            J = 4, # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200,
                            verbose = FALSE
  )
  pp
  expect_equal( nrow(pp), 2 * 3 * 2)


  grid <- pump_power_grid( d_m="d3.2_m3ff2rc",
                           MTP = "BF",
                           MDES = 0.1,
                           M = 3,
                           J = 3, # number of schools/block
                           K = 21, # number RA blocks
                           nbar = 258,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, numCovar.2 = 3,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = seq( 0, 0.3, 0.05 ),
                           ICC.3 = seq( 0, 0.45, 0.15 ),
                           rho = 0.4,
                           tnum = 100 )
  a <- length( unique( grid$ICC.2 ) )
  b <- length( unique( grid$ICC.3 ) )
  expect_equal( nrow( grid ), a * b * 2 )
  expect_true( "MDES" %in% names(grid) )


  grid2 <- pump_power_grid( d_m="d3.2_m3ff2rc", drop.unique.columns = FALSE,
                            MTP = "BF",
                            MDES = 0.1,
                            M = 3,
                            J = c( 3, 4 ), # number of schools/block
                            K = 21, # number RA blocks
                            nbar = 258,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.05, # significance level
                            numCovar.1 = 5, numCovar.2 = 3,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0,
                            ICC.3 = 0,
                            rho = 0.4,
                            tnum = 100 )
  grid2
  expect_true( "rho" %in% names( grid2 ) )

  grid3 <- pump_power_grid( d_m=c( "d3.2_m3ff2rc", "d3.2_m3rr2rc" ), drop.unique.columns = FALSE,
                            MTP = "BF",
                            MDES = 0.1,
                            M = 3,
                            J = 3, # number of schools/block
                            K = 21, # number RA blocks
                            nbar = 258,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.05, # significance level
                            numCovar.1 = 5, numCovar.2 = 3,
                            R2.1 = c( 0.1, 0.2 ), R2.2 = 0.7,
                            ICC.2 = 0,
                            ICC.3 = 0,
                            rho = c( 0, 0.4 ),
                            tnum = 100 )
  grid3
  expect_true( length( unique( grid3$d_m ) ) == 2 )



})


test_that("pump_mdes_grid works", {
    
  skip_on_cran()

  pp <- pump_mdes_grid(    d_m = "d3.2_m3ff2rc",
                           MTP = "HO",
                           target.power = c( 0.50, 0.80 ),
                           power.definition = "D1indiv",
                           tol = 0.05,
                           M = 5,
                           J = c( 3, 9), # number of schools/block
                           K = 7, # number RA blocks
                           nbar = 58,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.15, # significance level
                           numCovar.1 = 1, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = 0.05, ICC.3 = 0.9,
                           rho = 0.4, # how correlated outcomes are
                           verbose = FALSE, tnum = 1000,
  )
  pp
  expect_equal( nrow(pp), 2 * 2)

})


test_that("pump_sample_grid works", {
    
  skip_on_cran()

  pp <- pump_sample_grid(    d_m = "d3.2_m3ff2rc",
                             typesample = "J",
                             MTP = "HO",
                             MDES = 0.10,
                             target.power = c( 0.50, 0.80 ),
                             power.definition = "min1",
                             tol = 0.03,
                             M = 5,
                             K = 7, # number RA blocks
                             nbar = 58,
                             Tbar = 0.50, # prop Tx
                             alpha = 0.15, # significance level
                             numCovar.1 = 1, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = 0.25, ICC.3 = 0.25,
                             rho = 0.4, # how correlated outcomes are
                             verbose = FALSE, tnum = 400

  )
  pp
  expect_equal( nrow(pp), 2 )

})


test_that( "grid allows multiple MTP and power definitions", {

  skip_on_cran()
  pp <- pump_mdes_grid(    d_m = "d3.2_m3ff2rc",
                           MTP = c( "BF", "HO" ),
                           target.power = 0.5,
                           power.definition = c( "min1", "D1indiv" ),
                           tol = 0.05,
                           M = 5,
                           J = 5,
                           K = 7, # number RA blocks
                           nbar = 58,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.15, # significance level
                           numCovar.1 = 1, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = 0.05, ICC.3 = 0.9,
                           rho = 0.4, # how correlated outcomes are
                           verbose = FALSE, tnum = 500,
  )
  pp
  expect_equal( nrow( pp ), 4 )



  pp <- pump_sample_grid(    d_m = "d3.2_m3ff2rc",
                             MTP = c( "BF", "HO" ),
                             target.power = 0.9,
                             typesample = "J",
                             MDES = 0.05,
                             power.definition = c( "min1", "D1indiv" ),
                             tol = 0.05,
                             M = 5,
                             nbar = 10,
                             K = 12, # number RA blocks
                             Tbar = 0.50, # prop Tx
                             alpha = 0.15, # significance level
                             numCovar.1 = 1, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = 0.05, ICC.3 = 0.9,
                             rho = 0.4, # how correlated outcomes are
                             verbose = FALSE, tnum = 500,
  )
  pp
  expect_true( nrow( pp ) == 4 )
})




test_that( "grid works for long tables", {
  skip_on_cran()    
  pp <- pump_power_grid(    d_m = "d3.2_m3ff2rc",
                            MTP = c( "HO", "BF" ),
                            MDES = 0.10,
                            J = c( 4, 8 ),
                            M = 5,
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.25, ICC.3 = 0.25,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200, verbose = FALSE,
                            long.table = TRUE )

  expect_true( sum( is.na( pp$None ) ) > 0 )
  expect_true( ncol( pp ) == 7 )
})


test_that( "grid breaks with invalid inputs", {
  skip_on_cran()
  expect_error(pp <- pump_sample_grid(
    d_m = "d2.2_m2rc",
    MTP = c("HO"),
    typesample = c("J"),
    MDES = 0.125,
    J = 10,
    M = 5,
    nbar = 50,
    target.power = 0.8,
    power.definition = "indiv.mean",
    alpha = 0.5,
    Tbar = 0.8,
    numCovar.1 = 2,
    rho = 0.2))

  expect_error(pp <- pump_mdes_grid(
    d_m = "d2.2_m2rc",
    MTP = c("HO"),
    MDES = c(0.125),
    J = 10,
    M = 5,
    nbar = 50,
    target.power = 0.8,
    power.definition = "indiv.mean",
    alpha = 0.5,
    Tbar = 0.8,
    numCovar.1 = 2,
    rho = 0.2))
  
  expect_error(pp <- pump_mdes_grid(
    d_m = "d2.2_m2rc",
    MTP = c("HO"),
    J = c(10, 30),
    M = 5,
    nbar = c( 50, 100, 50 ),
    target.power = 0.8,
    power.definition = "indiv.mean",
    alpha = 0.5,
    Tbar = 0.8,
    numCovar.1 = 2, tnum = 100,
    rho = 0.2, tol = 0.5))
  
  expect_error(pp <- pump_power_grid(
      d_m = "d2.2_m2rc",
      MTP = c("HO"),
      MDES = 0.125,
      J = 10,
      M = 5,
      nbar = 50,
      alpha = 0.5,
      Tbar = 0.8,
      numCovar.1 = 2,
      rho = c(0.2, 0.2)))
  
})



test_that( "0 MDES does something reasonable", {
    
 
    pp1 = pump_power( d_m = "d1.1_m1c",
                      MDES = 0.3, M = 1, Tbar = 0.5,
                      nbar = 100 )
    pp1
    
    pp2 = update_grid( pp1, MDES = seq( 0, 1, by=0.1 ) )
    pp2

    expect_equal( nrow( pp2 ), 11 )
    expect_equal( pp2$D1indiv[[1]], 0.05 )
    
} )
