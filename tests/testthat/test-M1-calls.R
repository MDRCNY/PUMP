# library( PUMP )
# library( testthat )



test_that("M=1 calls on pump_power work", {
    
    set.seed( 2424424 )
    
    
    expect_warning( p1 <- pump_power( d_m = "d2.1_m2fc", 
                     MDES = 0.2,
                     Tbar = 0.5,
                     MTP = "BH",
                     nbar = 50, J = 20 ) )
    p1

    p2 <- pump_power( d_m = "d2.1_m2fc", 
                                      MDES = 0.2,
                                      Tbar = 0.5,
                                      nbar = 50, J = 20 )
    
    
    p3 <- pump_power( d_m = "d2.1_m2fc", 
                      M = 2,
                      MTP = "BF",
                      MDES = 0.2,
                      Tbar = 0.5,
                      nbar = 50, J = 20, rho = 0, tnum = 10000 )
    
    expect_true( p3$D1indiv[[1]] - p2$D1indiv < 0.01 )
    
    expect_equal( p1, p2 )
    
    # Should we be doing exact calculations when M=1?  Yes... but we
    # don't currently?
    expect_true( PUMP:::exact_calc(p2) )
    
} )



test_that("M=1 calls work", {
    
    set.seed( 2424424 )
    
    
    p1 = pump_mdes( d_m = "d2.1_m2fc", 
                    M = 1,
                    target.power = 0.80,
                    Tbar = 0.5,
                    nbar = 50, J = 20 )
    expect_true( !is.null( p1$`D1indiv power` ) )
    
    
    # Should this give an error?  Perhaps not?
    p1 = pump_mdes( d_m = "d2.1_m2fc", 
                    M = 1,
                    target.power = 0.80,
                    power.definition = "BH",
                    Tbar = 0.5,
                    nbar = 50, J = 20 )
    p1
    
    p1 = pump_mdes( d_m = "d2.1_m2fc", 
                    M = 1,
                    target.power = 0.80,
                    power.definition = "BO",
                    Tbar = 0.5,
                    nbar = 50, J = 20 )
    p1
    
    
    expect_error( p1 = pump_mdes( d_m = "d2.1_m2fc", 
                                  M = 2,
                                  target.power = 0.80,
                                  power.definition = "BH",
                                  Tbar = 0.5,
                                  nbar = 50, J = 20 ) )
    
    
} )
