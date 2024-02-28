# library( PUMP )
# library( testthat )
library(PowerUpR)

# NOTE: This should be even closer, no?  Why not perfect match?
default.tol <- 0.002




test_that( "very simple RCT tests", {
    sink("sink.txt")
    
    set.seed( 449595 )
    
    # Very large, low df, no covariates
    powerup.power <- PowerUpR::power.ira(
        es = 0.15,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        n = 1000
    )
    
    pump.power <- pump_power(
        d_m = "d1.1_m1c",
        MTP = 'None',
        nbar = 1000,
        M = 1,
        MDES = 0.15,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    )
    pump.power
    
    # So good so far, with high df the normal approximation holds
    expect_equal( powerup.power$df, pump.power$df1[[1]] )
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    
    
    # Very small, low df, no covariates
    powerup.power <- PowerUpR::power.ira(
        es = 1.25,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        n = 10
    )
    
    pump.power <- pump_power(
        d_m = "d1.1_m1c",
        MTP = 'None',
        nbar = 10,
        M = 1,
        MDES = 1.25,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    )
    
    # Via simulation
    pump.power.sim <- pump_power(
        d_m = "d1.1_m1c",
        MTP = 'None',
        nbar = 10,
        M = 1,
        MDES = 1.25,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE, 
        exact.where.possible = FALSE
    )
    
    pump.power
    pump.power.sim
    
    # Check df all line up
    expect_equal( powerup.power$df, pump.power$df1[[1]] )
    expect_equal( pump.power.sim$df1[[1]], pump.power$df1[[1]] )
    
    # Check power all line up
    expect_equal( pump.power.sim$D1indiv[1], pump.power$D1indiv[1], tol = default.tol )
    # TODO: Why this not align?
    #expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    
    #expect_equal(pump.power$D1indiv[1], pump.power.sim$D1indiv[1], tol = default.tol)
    
    
    
    # Very small, low df (test #2)
    powerup.power <- PowerUpR::power.ira(
        es = 1.25,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g1 = 5,
        r21 = 0.5,
        n = 10
    )
    
    pump.power <- pump_power(
        d_m = "d1.1_m1c",
        MTP = 'None',
        nbar = 10,
        M = 1,
        MDES = 1.25,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 5,
        R2.1 = 0.5, rho = 0
    )
    
    pump.power.sim <- pump_power(
        d_m = "d1.1_m1c",
        MTP = 'None',
        nbar = 10,
        M = 1,
        MDES = 1.25,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 5,
        R2.1 = 0.5, rho = 0, exact.where.possible = FALSE, tnum = 100000
    )
    
    powerup.power$power
    pump.power
    pump.power.sim
    
    # Check df all line up
    expect_equal( powerup.power$df, pump.power$df1[[1]] )
    expect_equal( pump.power.sim$df1[[1]], pump.power$df1[[1]] )
    
    # Check power all line up
    expect_equal( pump.power.sim$D1indiv[1], pump.power$D1indiv[1], tol = default.tol )
    # TODO: Why this not align?
    #expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    
    #expect_equal(pump.power$D1indiv[1], pump.power.sim$D1indiv[1], tol = default.tol)
    
    
    
    sink()
    
})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d2.1_m2fc -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


test_that( "testing of d2.1_m2fc", {
    
    sink("sink.txt")
    
    powerup.power <- PowerUpR::power.bira2c1(
        es = 0.125,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g1 = 1,
        r21 = 0.1,
        n = 50,
        J = 30
    )
    
    pump.power <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = 'None',
        nbar = 50,
        J = 30,
        M = 1,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1,
        R2.1 = 0.1, ICC.2 = 0, rho = 0,
        tnum = 100000
    )
    
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    powerup.power <- PowerUpR::power.bira2c1(
        es = 0.125,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g1 = 1,
        r21 = 0.1,
        n = 50,
        J = 30
    )
    
    pump.power <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = 'None',
        nbar = 50,
        J = 30,
        M = 1,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1,
        R2.1 = 0.1, ICC.2 = 0, rho = 0,
        tnum = 100000
    )
    
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    
    
    
    
    sink()
    file.remove("sink.txt")
    
})


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d2.1_m2fr -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fr one-tailed", {
    
    sink("sink.txt")
    
    set.seed(23423432)
    
    powerup.power <- expect_warning(PowerUpR::power.bira2r1(
        es = 0.125,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        omega2 = 0.1,
        r21 = 0.3,
        r2t2 = 0,
        n = 50,
        J = 30
    ))
    
    pump.power <- pump_power(
        d_m = "d2.1_m2fr",
        MTP = 'None',
        nbar = 50,
        J = 30,
        M = 1,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1,
        R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
        omega.2 = 0.1,
        tnum = 100000
    )
    
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    powerup.mdes <- expect_warning(PowerUpR::mdes.bira2r1(
        power = 0.8,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        omega2 = 0.1,
        r21 = 0.3,
        r2t2 = 0,
        n = 50,
        J = 30
    ))
    
    pump.mdes <- pump_mdes(
        target.power = 0.8,
        power.definition = 'D1indiv',
        d_m = "d2.1_m2fr",
        MTP = 'None',
        nbar = 50,
        J = 30,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1,
        R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
        omega.2 = 0.1
    )
    
    expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)
    
    
    powerup.ss <- expect_warning(PowerUpR::mrss.bira2r1(
        power = 0.8,
        es = 0.125,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        omega2 = 0.1,
        r21 = 0.3,
        r2t2 = 0,
        n = 50
    ))
    
    pump.ss <- pump_sample(
        target.power = 0.8,
        power.definition = 'D1indiv',
        typesample = 'J',
        d_m = "d2.1_m2fr",
        MTP = 'None',
        MDES = 0.125,
        nbar = 50,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1,
        R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
        omega.2 = 0.1
    )
    
    expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)
    
    sink()
    file.remove("sink.txt")
    
})


test_that("testing of d2.1_m2fr two-tailed", {
    
    sink("sink.txt")
    
    powerup.power <- expect_warning(PowerUpR::power.bira2r1(
        es = 0.125,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        omega2 = 0.1,
        r21 = 0.3,
        r2t2 = 0,
        n = 50,
        J = 30
    ))
    
    pump.power <- pump_power(
        d_m = "d2.1_m2fr",
        MTP = 'None',
        nbar = 50,
        J = 30,
        M = 1,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1,
        R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
        omega.2 = 0.1,
        tnum = 100000
    )
    
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    powerup.mdes <- expect_warning(PowerUpR::mdes.bira2r1(
        power = 0.8,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        omega2 = 0.1,
        r21 = 0.3,
        r2t2 = 0,
        n = 50,
        J = 30
    ))
    
    pump.mdes <- pump_mdes(
        target.power = 0.8,
        power.definition = 'D1indiv',
        d_m = "d2.1_m2fr",
        MTP = 'None',
        nbar = 50,
        J = 30,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1,
        R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
        omega.2 = 0.1
    )
    
    expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)
    
    
    powerup.ss <- expect_warning(PowerUpR::mrss.bira2r1(
        power = 0.8,
        es = 0.125,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        omega2 = 0.1,
        r21 = 0.3,
        r2t2 = 0,
        n = 50
    ))
    
    pump.ss <- pump_sample(
        target.power = 0.8,
        power.definition = 'D1indiv',
        typesample = 'J',
        d_m = "d2.1_m2fr",
        MTP = 'None',
        MDES = 0.125,
        nbar = 50,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1,
        R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
        omega.2 = 0.1
    )
    
    expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)
    
    sink()
    file.remove("sink.txt")
    
})


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d3.2_m3ff2rc -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d3.2_m3ff2rc one-tailed", {
    
    sink("sink.txt")
    
    powerup.power <- PowerUpR::power.bcra3f2(
        es = 0.125,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        r21 = 0.3,
        r22 = 0.3,
        n = 50,
        J = 30,
        K = 10
    )
    
    pump.power <- pump_power(
        d_m = "d3.2_m3ff2rc",
        MTP = 'None',
        nbar = 50,
        J = 30,
        K = 10,
        M = 1,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
        tnum = 100000
    )
    
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    powerup.mdes <- PowerUpR::mdes.bcra3f2(
        power = 0.8,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        r21 = 0.3,
        r22 = 0.3,
        n = 50,
        J = 30,
        K = 10
    )
    
    pump.mdes <- pump_mdes(
        target.power = 0.8,
        power.definition = 'D1indiv',
        d_m = "d3.2_m3ff2rc",
        MTP = 'None',
        nbar = 50,
        J = 30,
        K = 10,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
    )
    
    expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)
    
    
    powerup.ss <- PowerUpR::mrss.bcra3f2(
        power = 0.8,
        es = 0.125,
        alpha = 0.05,
        two.tailed = FALSE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        r21 = 0.3,
        r22 = 0.3,
        n = 50,
        J = 30,
        K = 10
    )
    
    pump.ss <- pump_sample(
        target.power = 0.8,
        power.definition = 'D1indiv',
        typesample = 'J',
        d_m = "d3.2_m3ff2rc",
        MTP = 'None',
        MDES = 0.125,
        nbar = 50,
        K = 10,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
        tol = 0.005
    )
    
    print(pump.ss$Sample.size)
    print(powerup.ss$J)
    expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)
    
    sink()
    file.remove("sink.txt")
    
})

test_that("testing of d3.2_m3ff2rc two-tailed", {
    
    sink("sink.txt")
    
    powerup.power <- PowerUpR::power.bcra3f2(
        es = 0.125,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        r21 = 0.3,
        r22 = 0.3,
        n = 50,
        J = 30,
        K = 10
    )
    
    pump.power <- pump_power(
        d_m = "d3.2_m3ff2rc",
        MTP = 'None',
        nbar = 50,
        J = 30,
        K = 10,
        M = 1,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
        tnum = 100000
    )
    
    expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
    
    powerup.mdes <- PowerUpR::mdes.bcra3f2(
        power = 0.8,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        r21 = 0.3,
        r22 = 0.3,
        n = 50,
        J = 30,
        K = 10
    )
    
    pump.mdes <- pump_mdes(
        target.power = 0.8,
        power.definition = 'D1indiv',
        d_m = "d3.2_m3ff2rc",
        MTP = 'None',
        nbar = 50,
        J = 30,
        K = 10,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
    )
    
    expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)
    
    
    powerup.ss <- PowerUpR::mrss.bcra3f2(
        power = 0.8,
        es = 0.125,
        alpha = 0.05,
        two.tailed = TRUE,
        p = 0.5,
        g2 = 1,
        rho2 = 0.05,
        r21 = 0.3,
        r22 = 0.3,
        n = 50,
        J = 30,
        K = 10
    )
    
    pump.ss <- pump_sample(
        target.power = 0.8,
        power.definition = 'D1indiv',
        typesample = 'J',
        d_m = "d3.2_m3ff2rc",
        MTP = 'None',
        MDES = 0.125,
        nbar = 50,
        K = 10,
        M = 1,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
        tol = 0.005
    )
    
    expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)
    
    sink()
    file.remove("sink.txt")
    
})

