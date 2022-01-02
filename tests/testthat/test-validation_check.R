# library( PUMP )
# library( testthat )


test_that("validation works at least vaguely", {

    params.list <- list(
        MTP = "HO",
        M = 4,
        J = 44,
        nbar = 1000,
        MDES = rep(0.40, 4),
        Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
        rho = 0.2,  tnum=1000
    )
    
    chk <- validate_inputs(d_m = "d2.2_m2rc", params.list = params.list)
    expect_true(!is.null(chk))
    
    params.list <- list(
        MTP = "HO",
        M = 4,
        J = -3,
        nbar = 1000,
        MDES = rep(0.40, 4),
        Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
        rho = 0.2,  tnum=1000
    )
    
    expect_error(chk <- validate_inputs(d_m = "d2.2_m2rc", params.list = params.list))


    params.list <- list(
        MTP = "HO",
        M = 4,
        J = 3,
        nbar = 1000,
        MDES = 0.4,
        numZero = 2,
        Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
        R2.3 = 0.5, omega.3 = 0.4,
        rho = 0.2,  tnum=1000
    )
    
    expect_warning(expect_message(chk <- validate_inputs(d_m = "d2.2_m2rc", params.list = params.list)))
    
    
    params.list <- list(
        MTP = "HO",
        M = 4,
        J = -4,
        K = 4,
        nbar = 1000,
        MDES = 0.4,
        Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
        R2.3 = 0.5, omega.3 = 0.4,
        rho = 0.2,  tnum=1000, verbose = TRUE
    )
    
    expect_error(chk <- validate_inputs(d_m = "d2.2_m2rc", params.list = params.list))


    params.list <- list(
        MTP = "HO",
        call = "sample",
        J = 10,
        nbar = 200,
        M = 20,
        MDES = rep(0.05, 20),
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
        rho = 0.95 
    )
    
    expect_warning(chk <- validate_inputs(d_m = "d2.1_m2fc", params.list = params.list))


    params.list <- list(
        MTP = "HO",
        call = "mdes",
        J = 10,
        nbar = 200,
        M = 20,
        MDES = rep(0.05, 20),
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
        rho = 0.95 
    )
    
    expect_warning(chk <- validate_inputs(d_m = "d2.1_m2fc", params.list = params.list))


    params.list <- list(
        MTP = "HO",
        call = "mdes",
        J = 10,
        nbar = 200,
        M = 20,
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = -0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
        rho = 0.95
    )
    
    expect_error(chk <- validate_inputs(d_m = "d2.1_m2fc", params.list = params.list))
    

    params.list <- list(
        MTP = c( "HO", "BF" ),
        call = "mdes",
        J = 10,
        nbar = 200,
        M = 20,
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.9, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
        rho = 0.95
    )
    
    expect_error(chk <- validate_inputs(d_m = "d2.1_m2fc", params.list = params.list))

} )
