# library( PUMP )
# library( testthat )

skip_on_cran()

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
  
  # Basic single vary call
  grid <- pump_power_grid( d_m = "d2.1_m2fr",
                           MTP = c( "BF", "BH" ),
                           MDES = c( 0.1, 0.2, 0.3 ),
                           M = 3,
                           J = 10, 
                           nbar = 258,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, 
                           R2.1 = 0.1,
                           ICC.2 = 0.15,
                           rho = 0.9,
                           tnum = 100,
                           long.table = FALSE )
  expect_true( length( attr( grid, "var_names") ) == 1 )
  expect_true( attr( grid, "var_names") == "MDES" )
  
  grid.plot <- plot(grid, power.definition = 'min1' )
  expect_true(!is.null(grid.plot))
  
  grid.plot <- plot(grid, power.definition = 'D2indiv' )
  expect_true(!is.null(grid.plot))
  
  grid.plot <- plot(grid)
  expect_true(!is.null(grid.plot))
  
  
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
  
  expect_true( length( attr( grid, "var_names") ) == 2 )
  
  grid.plot <- plot(grid, var.vary = 'ICC.3')
  expect_true(!is.null(grid.plot))
  grid.plot <- plot(grid, power.definition = 'indiv.mean', var.vary = 'ICC.2')
  expect_true(!is.null(grid.plot))
  
  # works for other definitions of power
  grid.plot <- plot(grid, power.definition = 'D1indiv', var.vary = 'ICC.3')
  expect_true(!is.null(grid.plot))
  
  
  # Check single outcome case, and longtable is false case
  grid <- pump_power_grid( d_m = "d3.2_m3ff2rc",
                           MTP = c("None"),
                           MDES = 0.1,
                           M = 1,
                           J = 3, # number of schools/block
                           K = 21, # number RA blocks
                           nbar = c( 120, 258 ),
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, numCovar.2 = 3,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = c( 0, 0.15, 0.3 ),
                           ICC.3 = 0.2,
                           rho = 0.4,
                           tnum = 100,
                           long.table = FALSE)

  grid.plot <- plot(grid, power.definition = 'D1indiv', var.vary = 'ICC.2')
  expect_true(!is.null(grid.plot))
  grid.plot <- plot(grid, var.vary = 'ICC.2')
  expect_true(!is.null(grid.plot))
  
})



test_that("Grid plot works for MDES", {
    
    grid <- pump_mdes_grid(  d_m = "d3.2_m3ff2rc",
                             MTP = c("BF", "BH"),
                             target.power = 0.8,
                             power.definition = "min1",
                             M = 3,
                             J = 3, # number of schools/block
                             K = 21, # number RA blocks
                             nbar = 258,
                             Tbar = 0.50, # prop Tx
                             alpha = 0.05, # significance level
                             numCovar.1 = 5, numCovar.2 = 3,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = 0.4,
                             ICC.3 = seq( 0, 0.45, 0.15 ),
                             rho = 0.4, tnum = 100, tol = 0.45 )
    
    expect_true( length( attr( grid, "var_names") ) == 2 )
    
    gg <- handle_power_definition(grid, "min1", "MDES", "ICC.3", FALSE )
    expect_true( gg$powerType == "1-minimum" )
    expect_true( !gg$multiPower)
    expect_true( is.null( gg$title ) )
    
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
                             K = c( 10, 20, 30 ), # number RA blocks
                             nbar = 258,
                             Tbar = 0.50, # prop Tx
                             alpha = 0.05, # significance level
                             numCovar.1 = 5, numCovar.2 = 3,
                             R2.1 = 0.1, R2.2 = 0.7,
                             ICC.2 = 0.15,
                             ICC.3 = 0.4,
                             rho = 0.4, tnum = 100, tol = 0.45 ))
    grid.plot <- plot(grid, power.definition = 'complete', var.vary = 'K')
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
        rho = seq(0.2, 0.8, 0.2), 
        tnum = 100, tol = 0.45 ))
    grid.plot <- plot(grid)
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
  
  
  grid <- pump_sample_grid(  d_m = "d2.2_m2rc",
                                             MTP = c( "HO", "BH" ),
                                             target.power = 0.8,
                                             power.definition = c( "min1", 'complete' ),
                                             typesample = 'J',
                                             MDES = 0.1,
                                             M = 3,
                                             nbar = 258,
                                             Tbar = 0.50, # prop Tx
                                             alpha = 0.05, # significance level
                                             numCovar.1 = 5, numCovar.2 = 3,
                                             R2.1 = 0.1, R2.2 = 0.7,
                                             ICC.2 = c( 0, 0.3 ),
                                             rho = 0.4, tnum = 100, tol = 0.45 )
  # grid
  grid.plot <- plot(grid, power.definition = 'complete' )
  expect_true(!is.null(grid.plot))
  
  # grid.plot (split plot)
  grid.plot <- plot(grid )
  expect_true(!is.null(grid.plot))

  
})


test_that( "power_curve works", {
  
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


test_that( "power curve plotting works", {
    
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
    
    sp <- search_path(nbar3)
    expect_true(!is.null(sp))
    p <- PUMP:::plot_power_search(nbar3)
    expect_true(!is.null(p))
    
    
    pts <- PUMP:::get_sample_tick_marks(desired_pts = sp$pt[sp$pt<1000], 
                                        breaks = 10, 
                                       include.points = TRUE, log = TRUE )
    expect_true( length( pts ) == 10 )
    pts <- PUMP:::get_sample_tick_marks(desired_pts = sp$pt[sp$pt<1000], breaks = 10, 
                                       include.points = TRUE, log = FALSE )
    expect_true( length( pts ) == 10 )
    
    
    expect_true(!is.null(p <- power_curve(nbar3)))  
    expect_true(!is.null(p <- plot_power_curve(nbar3)))  
    
    mdes <- expect_warning(pump_mdes(d_m = "d2.1_m2fc", MTP = 'HO',
                      power.definition = 'D1indiv', 
                      target.power = 0.4,
                      J = 10, nbar = 20, M = 3, 
                      Tbar = 0.5, alpha = 0.05,
                      numCovar.1 = 1, R2.1 = 0.1, 
                      ICC.2 = 0.05, rho = 0.2,
                      max.steps = 3))
    plot_power_curve(mdes)
})



