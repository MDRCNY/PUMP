


check_design <- function( d_m, K, J, nbar, dfstar,
                          numCovar.2 = 2, numCovar.3 = 4 ) {
    SE <- PUMP:::calc_SE(d_m = d_m, J = J, K = K, nbar = nbar, Tbar = 0.2, 
                         R2.1 = 0.3, R2.2 = 0.2, R2.3 = 0.4, ICC.2 = 0.15, ICC.3 = 0.25, 
                         omega.2 = 0.075, omega.3 = 0.275 )
    SE
    df <- PUMP:::calc_df(d_m = d_m, J = J, K = K, nbar = nbar,
                         numCovar.1 = 1, numCovar.2 = numCovar.2, numCovar.3 = numCovar.3 )
    expect_equal( df,  dfstar )
    
    nbarstar <- PUMP:::calc_nbar(d_m = d_m, MDES = 2.8*SE, J = J, K = K, Tbar = 0.2, 
                                 R2.1 = 0.3, R2.2 = 0.2, R2.3 = 0.4, ICC.2 = 0.15, ICC.3 = 0.25, 
                                 omega.2 = 0.075, omega.3 = 0.275 )
    expect_equal( nbar, nbarstar )
    
    Jstar <- PUMP:::calc_J(d_m = d_m, MDES = 2.8*SE, nbar = nbar, K = K, Tbar = 0.2, 
                           R2.1 = 0.3, R2.2 = 0.2, R2.3 = 0.4, ICC.2 = 0.15, ICC.3 = 0.25, 
                           omega.2 = 0.075, omega.3 = 0.275 )
    expect_equal( J, Jstar )
    
    
    Kstar <- PUMP:::calc_K(d_m = d_m, MDES = 2.8*SE, nbar = nbar, J = J, Tbar = 0.2, 
                           R2.1 = 0.3, R2.2 = 0.2, R2.3 = 0.4, ICC.2 = 0.15, ICC.3 = 0.25, 
                           omega.2 = 0.075, omega.3 = 0.275 )
    expect_equal( K, Kstar )
    
}


test_that("calc_SE and calc_J correspond", {
    
    K = 4
    J = 10
    nbar = 20
    
    check_design( d_m = "d3.1_m3ff2rr", K = K, J = J, nbar = nbar,
                  dfstar = K * J - K - 1)
    
    check_design( d_m = "d3.1_m3rr2rr", K = K, J = J, nbar = nbar,
                  dfstar = K - 1)
    
    
    numCovar.2 = 2
    check_design( d_m = "d3.2_m3ff2rc", K = K, J = J, nbar = nbar,
                  dfstar = K * (J-2) - numCovar.2 )
    
    check_design( d_m = "d3.2_m3rr2rc", K = K, J = J, nbar = nbar,
                  dfstar = K - 1 )
    
    numCovar.3 = 4
    check_design( d_m = "d3.3_m3rc2rc", K = 20, J = J, nbar = nbar,
                  dfstar = 20 - numCovar.3 - 2 )

    
})
