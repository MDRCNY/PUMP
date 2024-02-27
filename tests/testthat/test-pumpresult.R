# library( PUMP )
# library( testthat )

test_that( "pumpresult dimensions work", {

    set.seed( 101010 )
    ss <- pump_power(    d_m = "d2.1_m2fc",
                          MTP = c( "HO", "BH" ),
                          nbar = 200, J = 40,
                          M = 5,
                          MDES = 0.05,
                          Tbar = 0.50, alpha = 0.05,
                          numCovar.1 = 5,
                          R2.1 = 0.1,
                          ICC.2 = 0.05,
                          rho = 0, tnum = 100 )
    
    
    
    ss
    
    expect_equal( nrow(ss), 3 )

    expect_equal( ss[[1]], c( "None", "HO", "BH" ) )
    
    #expect_true( length(ss[2,] ) == 12 )
    expect_true( is.character(ss[[3,1]] ) )
    expect_true( all( is.na( ss[1,9:12] ) ) )
    
    ss
    ssL <- transpose_power_table(ss)
    ssL
    expect_true( is.pumpresult(ssL) )
    
    expect_equal( d_m( ssL ), "d2.1_m2fc" )

    
    ssLL <- transpose_power_table(ssL)   
    ssLL
    ss
    expect_equal( nrow( ssLL ), nrow( ss ) )
    expect_equal( colnames(ssLL), colnames(ss)[1:ncol(ssLL)] )
    expect_true( is.pumpresult(ssLL) )
    
    
    ssLt <- pump_power(    d_m = "d2.1_m2fc",
                           MTP = c( "HO", "BH" ),
                           nbar = 200, J = 40,
                           M = 5,
                           MDES = 0.05,
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 5,
                           R2.1 = 0.1,
                           ICC.2 = 0.05,
                           rho = 0, tnum = 100, long.table = TRUE )

    ssLtL <- transpose_power_table(ssLt)
    expect_true( is.pumpresult(ssLtL) )
    
    expect_equal( dim( ssLt ), dim( ssL ) )
    expect_equal( colnames( ssLt ), colnames( ssL ) )
    
})







test_that( "pumpresult for sample and mdes work", {
    
    set.seed( 10101033 )
    ss <- pump_sample(   d_m = "d2.1_m2fc",
                         MTP = "BH",
                         nbar = 200,
                         typesample = "J",
                         power.definition = "min1",
                         target.power = 0.80,
                         M = 5,
                         MDES = 0.05,
                         Tbar = 0.50, alpha = 0.05,
                         numCovar.1 = 5,
                         R2.1 = 0.1,
                         ICC.2 = 0.05,
                         rho = 0, start.tnum = 100, tnum = 100,
                         tol = 0.03 )
    
    
    dim(ss)
    
    expect_equal( dim(ss), c(1,4) )
    
    
})


