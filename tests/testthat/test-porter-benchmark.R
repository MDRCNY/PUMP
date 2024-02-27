


test_that( "testing against Porter 2017", {
    
    skip_on_cran()
    
    # NOTE: compare to Figure 1
    
    # all outcomes have effects
    # 3 outcomes
    # rho = 0.2
    set.seed(13434)
    pp1 <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = c( 'BF', 'HO', 'BH', 'WY-SS', 'WY-SD' ),
        nbar = 50,
        J = 20,
        M = 3,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, tnum = 1000, B = 1000,
        R2.1 = 0.5, ICC.2 = 0, rho = 0.2)
    
    expect_equal(0.65, pp1$D1indiv[2], tol = 0.01)
    expect_equal(0.72, pp1$D1indiv[3], tol = 0.01)
    expect_equal(0.77, pp1$D1indiv[4], tol = 0.01)
    expect_equal(0.68, pp1$D1indiv[5], tol = 0.01)
    expect_equal(0.75, pp1$D1indiv[6], tol = 0.01)
    
    # all outcomes have effects
    # 3 outcomes
    # rho = 0.8
    set.seed(13434)
    pp2 <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = c('BF', 'HO', 'BH', 'WY-SS', 'WY-SD'),
        nbar = 50,
        J = 20,
        M = 3,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, tnum = 1000, B = 1000,
        R2.1 = 0.5, ICC.2 = 0, rho = 0.8)
    
    expect_equal(0.65, pp2$D1indiv[2], tol = 0.01)
    expect_equal(0.71, pp2$D1indiv[3], tol = 0.01)
    expect_equal(0.76, pp2$D1indiv[4], tol = 0.01)
    expect_equal(0.71, pp2$D1indiv[5], tol = 0.01)
    expect_equal(0.76, pp2$D1indiv[6], tol = 0.01)
    
    
    # 2/3 outcomes have effects (1/3 are zero)
    # 9 outcomes
    # rho = 0.2
    set.seed(13434)
    pp3 <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = c('BF', 'HO', 'BH', 'WY-SS', 'WY-SD'),
        propZero = 1/3,
        nbar = 50,
        J = 20,
        M = 9,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, tnum = 1000, B = 1000,
        R2.1 = 0.5, ICC.2 = 0, rho = 0.2)
    
    expect_equal(0.48, pp3$D1indiv[2], tol = 0.01)
    expect_equal(0.56, pp3$D1indiv[3], tol = 0.01)
    expect_equal(0.69, pp3$D1indiv[4], tol = 0.01)
    expect_equal(0.52, pp3$D1indiv[5], tol = 0.01)
    expect_equal(0.55, pp3$D1indiv[6], tol = 0.01)
    
    
    # all outcomes have effects
    # 12 outcomes
    # rho = 0.8
    set.seed(13434)
    pp4 <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = c('BF', 'HO', 'BH', 'WY-SS', 'WY-SD'),
        nbar = 50,
        J = 20,
        M = 12,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, tnum = 1000, B = 1000,
        R2.1 = 0.5, ICC.2 = 0, rho = 0.8)
    
    expect_equal(0.48, pp4$D1indiv[2], tol = 0.01)
    expect_equal(0.57, pp4$D1indiv[3], tol = 0.01)
    expect_equal(0.72, pp4$D1indiv[4], tol = 0.01)
    expect_equal(0.60, pp4$D1indiv[5], tol = 0.01)
    expect_equal(0.68, pp4$D1indiv[6], tol = 0.01)
    
})
