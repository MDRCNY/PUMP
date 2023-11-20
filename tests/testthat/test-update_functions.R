# library( PUMP )
# library( testthat )

test_that( "update generally works", {

    set.seed( 101010 )
    ss <- pump_sample(    d_m = "d2.1_m2fc",
                          MTP = "HO",
                          typesample = "J",
                          nbar = 200,
                          power.definition = "min1",
                          M = 5,
                          MDES = 0.05, target.power = 0.8,
                          tol = 0.05,
                          Tbar = 0.50, alpha = 0.05,
                          numCovar.1 = 5,
                          R2.1 = 0.1,
                          ICC.2 = 0.15,
                          rho = 0,
                          final.tnum = 1000 )

    up <- update(ss, nbar = 10, tnum = 5000 )
    
    expect_true( up$`Sample.type` == ss$`Sample.type` )

    pm <- params( up )
    upm <- params( up )
    expect_true( length(pm) == length(upm) )
    expect_true( all( pm$MDES == upm$MDES ) )
    expect_true( ss$`Sample.size` < up$`Sample.size`)

    topow <- update( up, type = "power", J = 20 )
    expect_true( topow$min1[[2]] < up$`min1 power` )


    set.seed( 4440404 )
    tomdes <- update( topow, type="mdes",
                     target.power = 0.95,
                     power.definition = "min2", tol=0.02, tnum=500 )
    expect_true( tomdes$`Adjusted.MDES` > params( topow )$MDES[[1]] )

    bsamp <- update( tomdes, type="sample", typesample = "nbar",
                    MDES = 0.10 )
    expect_true( bsamp$`Sample.size` > params(tomdes)$nbar )

    tp2 <- update( topow, d_m = "d3.2_m3ff2rc", K = 10,
                  MDES = 0.02, ICC.3 = 0.4, R2.2 = 0.1, numCovar.2 = 1)
    tp2
    expect_true( d_m(tp2) == "d3.2_m3ff2rc" )
})



test_that( "unused parameters get preserved", {
    res <- pump_mdes( "d2.1_m2fr", nbar = 80, J = 23,
                       Tbar = 0.50,
                       R2.1 = 0.60, 
                       ICC.2 = 0.20, numCovar.1 = 5,
                       omega.2 = 0.4,
                       target.power = 0.80 )
    res2 = update( res, d_m = "d2.1_m2fc" )
    p = params(res2)
    expect_equal( p$omega.2, 0.4 )
    res3 = update( res2, d_m = "d2.1_m2fr" )
    p2 = params(res3)
    expect_equal( p2$omega.2, 0.4 )
    expect_equal( res, res3 )
    
    grd = update_grid( res, ICC.2 = c( 0.2, 0.4 ) )
    summary( grd )
    grd2 = update_grid( grd, d_m = "d2.1_m2fc" )
    expect_equal( params(grd2)$omega.2, 0.4 )
} )

    
test_that( "update works for parallel", {
    
    skip_on_cran()
    set.seed( 101010 )
    pp <- pump_power( d_m = "d3.2_m3ff2rc",
                      MTP = "WY-SD",
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
                      tnum = 400
    )
    
    set.seed( 101010 )
    up <- update(pp, parallel.WY.cores = 2)
    
    expect_equal( up$D1indiv[2], pp$D1indiv[2], tol = 0.1 )

})



test_that( "update_grid generally works", {
    
    skip_on_cran()

    set.seed( 101010 )
    ss <- pump_mdes(    d_m = "d2.1_m2fc",
                        MTP = "HO",
                        nbar = 200, J = 20,
                        power.definition = "complete",
                        M = 3, target.power = 0.5,
                        tol = 0.02,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1,
                        ICC.2 = 0.05,
                        rho = 0,
                        final.tnum = 1000 )
    ss

    gd <- update_grid( ss, J = c( 10, 20, 30 ) )
    gd
    expect_true( nrow( gd ) == 3 )

    s2 <- update( ss, type="power", MDES = 0.10 )
    gd2 <- update_grid( s2, R2.1 = c( 0.3, 0.5 ) )
    expect_true( nrow( gd2 ) == 4 )

    s3 <- update( s2, type="sample", typesample="J",
                 power.definition = "D2indiv",
                 target.power = 0.70, tol = 0.03 )
    gd3 <- update_grid( s3, power.definition = c( "min1", "min2" )  )
    expect_true( nrow(gd3) == 2 )
} )


test_that( "updating grids possible", {
    res2 <- pump_mdes( "d2.1_m2fc", nbar = 80, J = 23,
                       Tbar = 0.50,
                       R2.1 = 0.60, 
                       ICC.2 = 0.20, numCovar.1 = 5,
                       omega.2 = 0.4,
                       target.power = 0.80 )
    
    gd3 = update_grid( res2, omega.2 = c( 0.2, 0.4, 0.6 ) )
    expect_equal( nrow(gd3), 3 )
    
    # update grid when we have a grid?
    # Two ways--both methods should update.
    gd4 <- update_grid( gd3, omega.2 = c( 1, 2, 3, 4 ) )
    expect_equal( nrow(gd4), 4 )
    expect_equal( params(gd4), params(gd3) )
    
    gd5 = update_grid( gd3, ICC.2 = c( 0.1, 0.2, 0.3 ) )
    expect_equal( nrow( gd5 ), nrow( gd3 ) * 3 )
    p3 = params( gd3 )
    p5 = params( gd5 )
    p3$ICC.2 = NULL
    p5$ICC.2 = NULL
    expect_equal( p3, p5 )
    
    
    gd4b <- update( gd3, omega.2 = c( 1, 2, 3, 4 ) )
    expect_equal( gd4, gd4b )
} )




