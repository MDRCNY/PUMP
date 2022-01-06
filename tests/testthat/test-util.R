

test_that("get_power_names works", {

    pn = PUMP:::get_power_names(M=3)
    expect_true( length( pn ) == 7 )
    
    pn = PUMP:::get_power_names(M=1)
    pn
    expect_true( length( pn ) == 1 )
    
})
