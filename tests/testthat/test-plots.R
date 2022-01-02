# library( PUMP )
# library( testthat )

test_that("Single Scenario plot works", {

  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = c( "BF", "BH" ),
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 21, # number RA blocks
                    nbar = 258,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, # how correlated outcomes are
                    tnum = 200
  )
  
  ss.plot <- plot(pp)
  expect_true(!is.null(ss.plot))

  
  ppL <- update( pp, long.table=TRUE )
  
  ss.plot <- plot(ppL)
  #ss.plot
  expect_true(!is.null(ss.plot))
  
    
})



test_that("Grid plot works for power", {

  grid <- pump_power_grid( d_m = "d3.2_m3ff2rc",
                           MTP = c( "BF", "BH" ),
                           MDES = 0.1,
                           M = 3,
                           J = 3, # number of schools/block
                           K = 21, # number RA blocks
                           nbar = 258,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, numCovar.2 = 3,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = c( 0, 0.3 ),
                           ICC.3 = seq( 0, 0.45, 0.15 ),
                           rho = 0.4,
                           tnum = 100,
                           long.table = TRUE)
  
  grid.plot <- plot(grid, power.definition = 'min1', var.vary = 'ICC.3')
  expect_true(!is.null(grid.plot))
  # works for other definitions of power
  grid.plot <- plot(grid, power.definition = 'D1indiv', var.vary = 'ICC.3')
  expect_true(!is.null(grid.plot))
  
  grid <- pump_power_grid( d_m = "d3.2_m3ff2rc",
                           MTP = c("None"),
                           MDES = 0.1,
                           M = 1,
                           J = 3, # number of schools/block
                           K = 21, # number RA blocks
                           nbar = 258,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, numCovar.2 = 3,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = c( 0, 0.3 ),
                           ICC.3 = seq( 0, 0.45, 0.15 ),
                           rho = 0.4,
                           tnum = 100,
                           long.table = TRUE)
  
  grid.plot <- plot(grid, power.definition = 'D1indiv', var.vary = 'ICC.3')
  expect_true(!is.null(grid.plot))
})


test_that("Grid plot works for MDES", {
    
    grid <- pump_mdes_grid(  d_m = "d3.2_m3ff2rc",
                             MTP = c("BF", "BH"),
                             target.power = 0.8,
                             power.definition = 'min1',
                             M = 3,
                             J = 3, # number of schools/block
                             K = 21, # number RA blocks
                             nbar = 258,
                             Tbar = 0.50, # prop Tx
                             alpha = 0.05, # significance level
                             numCovar.1 = 5, numCovar.2 = 3,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = c( 0, 0.3 ),
                             ICC.3 = seq( 0, 0.45, 0.15 ),
                             rho = 0.4, tnum = 100, tol = 0.45 )
    
    grid.plot <- plot(grid, power.definition = 'min1', var.vary = 'ICC.3')
    expect_true(!is.null(grid.plot))

})




test_that("Grid plot works for SS", {
    
    grid <-expect_warning(pump_sample_grid(  d_m = "d3.2_m3ff2rc",
                             MTP = c( "HO", "BH" ),
                             target.power = 0.8,
                             power.definition = 'complete',
                             typesample = 'J',
                             MDES = 0.2,
                             M = 3,
                             K = 21, # number RA blocks
                             nbar = 258,
                             Tbar = 0.50, # prop Tx
                             alpha = 0.05, # significance level
                             numCovar.1 = 5, numCovar.2 = 3,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = c( 0, 0.3 ),
                             ICC.3 = seq( 0, 0.45, 0.15 ),
                             rho = 0.4, tnum = 100, tol = 0.45 ))
    grid.plot <- plot(grid, power.definition = 'complete', var.vary = 'ICC.3')
    expect_true(!is.null(grid.plot))
    
    grid <- expect_warning(pump_sample_grid(  
        d_m = "d3.2_m3ff2rc",
        MTP = c( "HO", "BH" ),
        target.power = 0.8,
        power.definition = 'D1indiv',
        typesample = 'J',
        MDES = 0.2,
        M = 3,
        K = 21, # number RA blocks
        nbar = 258,
        Tbar = 0.50, # prop Tx
        alpha = 0.05, # significance level
        numCovar.1 = 5, numCovar.2 = 3,
        R2.1 = 0.1, R2.2 = 0.7,
        ICC.2 = 0.2,
        ICC.3 = 0.2,
        rho = seq(0.2, 0.8, 0.1), 
        tnum = 100, tol = 0.45 ))
    grid.plot <- plot(grid, power.definition = 'D1indiv', var.vary = 'rho')
    expect_true(!is.null(grid.plot))
    
})


test_that("Two variable plot works for SS", {
  
  expect_warning( grid <- pump_sample_grid(  d_m = "d3.2_m3ff2rc",
                             MTP = c( "HO", "BH" ),
                             target.power = 0.8,
                             power.definition = 'complete',
                             typesample = 'K',
                             MDES = 0.1,
                             M = 3,
                             J = 21, # number RA blocks
                             nbar = 258,
                             Tbar = 0.50, # prop Tx
                             alpha = 0.05, # significance level
                             numCovar.1 = 5, numCovar.2 = 3,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = c( 0, 0.3 ),
                             ICC.3 = c( 0, 0.3, 0.6 ),
                             rho = 0.4, tnum = 100, tol = 0.45 ) )
  # grid
  grid.plot <- plot(grid, power.definition = 'complete', var.vary = 'ICC.3')
  # grid.plot
  
  expect_true(!is.null(grid.plot))
  
})


test_that( "power curve works", {
  
  set.seed( 101010 )
  up <- pump_sample(    d_m = "d2.1_m2fc",
                        MTP = "HO",
                        typesample = "J",
                        nbar = 10,
                        power.definition = "min1",
                        M = 5,
                        MDES = 0.05, target.power = 0.8,
                        tol = 0.05,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1,
                        ICC.2 = 0.05,
                        rho = 0,
                        final.tnum = 100 )
  
  pc <- power_curve( up, low = 5, high = 1000, tnum = 200, grid.size=20, all=TRUE )
  expect_true( is.data.frame(pc) )
  expect_true( nrow(pc) > 20 )
  
  pt <- plot_power_curve(pc)
  # pt
  expect_true( !is.null( pt ) )
  
  pc <- power_curve( up, low = 250, high = 1000, tnum = 100, grid.size=20, all=FALSE )
  expect_true( is.data.frame(pc) )
  expect_true( nrow(pc) == 20 )
  expect_true( all( pc$w == 100 ) )
})


test_that( "search functions work with non-convergence", {
    
    nbar3 <- expect_warning(pump_sample(
        d_m = "d3.3_m3rc2rc",
        power.definition = "D1indiv",
        target.power = 0.4,
        typesample = "nbar",
        MTP = "HO",
        K = 20,
        J = 40,
        M = 3,
        MDES = rep(0.25, 3),
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5
    ))
    
    expect_true(!is.null(search_path(nbar3)))
    expect_true(!is.null(expect_warning(p <- plot_power_search(nbar3))))
    expect_true(!is.null(p <- power_curve(nbar3)))  
    expect_true(!is.null(p <- plot_power_curve(nbar3)))  
})



