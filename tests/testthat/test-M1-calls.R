# library( PUMP )
# library( testthat )

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
