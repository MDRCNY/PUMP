
library( testthat )

test_that( "pumpgridresult dimensions work", {
    
    set.seed( 101010 )
    ss <- pump_power_grid(    d_m = "d2.1_m2fc",
                              MTP = c( "HO", "BH" ),
                              nbar = 200, J = 40,
                              M = 5,
                              numZero = 2,
                              MDES = c( 0.05, 0.10, 0.15 ),
                              Tbar = 0.50, alpha = 0.05,
                              numCovar.1 = c( 5, 10 ),
                              R2.1 = 0.1,
                              ICC.2 = 0.05,
                              rho = 0, tnum = 100 )
    
    atrs <- attributes(ss)
    names(atrs)
    expect_equal( attr(ss, "var_names" ), c( "MDES", "numCovar.1" ) )
    expect_true( is.pumpgridresult( ss ) )
    
    ss
    expect_equal( nrow(ss), 18 )

    expect_equal( length( ss[[1]] ), nrow(ss) )
    
    expect_true( length(ss[2,] ) == dim(ss)[[2]] )
    expect_true( is.character(ss[[3,1]] ) )
    
    #print(ss)
    capture_output( pp <- print( ss ) )
    
    capture_output( pd <- print_context(ss) )
    
    skip_on_cran()
    
    set.seed( 1010310 )
    ssL <- pump_power_grid(    d_m = "d2.1_m2fc",
                              MTP = c( "HO", "BH" ),
                              nbar = 200, J = 40,
                              M = 5,
                              MDES = c( 0.05, 0.10, 0.15 ),
                              Tbar = 0.50, alpha = 0.05,
                              numCovar.1 = 5,
                              R2.1 = 0.1,
                              ICC.2 = 0.05,
                              rho = 0, tnum = 100, long.table = TRUE )
    
    ssL
    expect_true( is.pumpgridresult( ssL ) )
    
    class(ssL)
    #expect_equal( dim(ssL), c(33, 6) )
    
    
    #ssL2 = transpose_power_table( ss )
    
    
    ssLW <- pump_mdes_grid(    d_m = "d2.1_m2fc",
                               nbar = 200, J = 40,
                               target.power = 0.80, 
                               M = 1,
                               Tbar = 0.50, alpha = 0.05,
                               numCovar.1 = 5,
                               R2.1 = c( 0.1, 0.4 ),
                               ICC.2 = 0.05,
                               rho = 0, tnum = 500, start.tnum = 100, 
                               tol = 0.5,
                               drop.unique.columns = FALSE)
    
    ssLW
    output <- capture.output(summary(ssLW))
    expect_true(!any(grepl("outcomes", output)))
    expect_equal( nrow(ssLW), 2)
    
    ssLW <- pump_mdes_grid(    d_m = "d2.1_m2fc",
                               MTP = c( "HO", "BH" ),
                               nbar = 200, J = 40,
                               target.power = c( 0.60, 0.70 ), 
                               power.definition="D2indiv",
                               M = 5,
                               Tbar = 0.50, alpha = 0.05,
                               numCovar.1 = 5,
                               R2.1 = c( 0.1, 0.4 ),
                               ICC.2 = 0.05,
                               rho = 0, start.tnum = 100,
                               tnum = 100, tol = 0.45,
                               drop.unique.columns = TRUE)
    
    ssLW
    
    set.seed( 40440 )
    ssLW <- pump_sample_grid(  d_m = "d3.2_m3rr2rc",
                               MTP = c( "HO", "BH" ),
                               nbar = 200, J = 40,
                               target.power = c( 0.60, 0.70 ), 
                               typesample = "K",
                               MDES = 0.1,
                               power.definition="D2indiv",
                               M = 5,
                               Tbar = 0.50, alpha = 0.05,
                               numCovar.1 = 5,
                               R2.1 = c( 0.1, 0.4 ),
                               ICC.2 = 0.05,
                               rho = 0, tnum = 100, tol = 0.45,
                               drop.unique.columns = TRUE)
    
    ssLW
    
    
    expect_true( is.pumpgridresult( ssLW ) )
    
    expect_output( ss <- summary( ssLW ) )
    expect_true( !is.null( ss ) )

    
})




