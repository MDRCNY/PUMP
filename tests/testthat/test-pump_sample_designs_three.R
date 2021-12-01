# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ----- three level models ------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


source( "testing_code.R" )


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# --------    d3.1_m3rr2rr    --------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d3.1_m3rr2rr", {

    if ( FALSE ) {

        set.seed( 524235326 )
        pp1 <- pump_power(
            design = "d3.1_m3rr2rr",
            MTP = 'Holm',
            nbar = 50,
            K = 15,
            J = 30,
            M = 3,
            MDES = rep(0.125, 3),
            Tbar = 0.5, alpha = 0.05,
            numCovar.1 = 1, numCovar.2 = 1,
            R2.1 = 0.1, R2.2 = 0.1,
            ICC.2 = 0.2, ICC.3 = 0.2,
            omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
            tnum = 100000)
        pp1
        pp1_power = pp_power
    }

    # Save the target power from above.  No need to rerun for testing code
    pp1_power = 0.79038

    vals = test_sample_triad( target_power = pp1_power,
                              seed = 524235326,
                              nbar = 50, J = 30, K = 15,
                              design = "d3.1_m3rr2rr",
                              MTP = 'Holm',
                              power.definition = 'D1indiv',
                              M = 3,
                              MDES = 0.125,
                              Tbar = 0.5, alpha = 0.05,
                              numCovar.1 = 1, numCovar.2 = 1,
                              R2.1 = 0.1, R2.2 = 0.1,
                              ICC.2 = 0.2, ICC.3 = 0.2,
                              omega.2 = 0.1, omega.3 = 0.1, rho = 0.5 )

    expect_equal(50, vals$nbar, tol=0.20)
    expect_equal(30, vals$J, tol=0.10)
    expect_equal(15, vals$K, tol=0.10)

    expect_equal( warning_pattern(vals), c(TRUE, TRUE,FALSE) )


    # if we go below the true value, we get the wrong number since it is so flat
    set.seed( 524235325 )

    expect_warning(nbar7 <- pump_sample(
        design = "d3.1_m3rr2rr",
        typesample = 'nbar',
        MTP = 'Holm',
        target.power = pp1_power,
        power.definition = 'D1indiv',
        K = 15,
        J = 30,
        M = 3,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.1,
        ICC.2 = 0.2, ICC.3 = 0.2,
        omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
        max_sample_size_nbar = 40 ))
    #  plot_power_search(nbar7)
    expect_true(nbar7$`Sample.size` < 45 )
})



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------ d3.2_m3ff2rc ------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


test_that("testing of d3.2_m3ff2rc", {
    set.seed( 245444 )

    if ( FALSE ) {
        pp1 <- pump_power(
            design = "d3.2_m3ff2rc",
            MTP = 'Holm',
            nbar = 50,
            J = 30,
            K = 10,
            M = 5,
            MDES = 0.125,
            Tbar = 0.5, alpha = 0.05,
            numCovar.1 = 1, numCovar.2 = 1,
            R2.1 = 0.1, R2.2 = 0.1,
            ICC.2 = 0.2, ICC.3 = 0.2,
            omega.2 = 0, omega.3 = 0.1, rho = 0.5, tnum = 100000)
        pp1
        pp1$min2[2]
    }
    pp_power = 0.64632

    # In this case it looks
    vals <- test_sample_triad( pp_power, nbar=50, J=30, K=10, 4224422,
                               design = "d3.2_m3ff2rc",
                               MTP = 'Holm',
                               power.definition = 'min2',
                               M = 5,
                               MDES = 0.125,
                               Tbar = 0.5, alpha = 0.05,
                               numCovar.1 = 1, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.1,
                               ICC.2 = 0.2, ICC.3 = 0.2,
                               omega.2 = 0, omega.3 = 0.1, rho = 0.5 )
    vals[1:3]

    # It is super flat for nbar!!!  Big tolerance.
    expect_equal(50, vals$nbar, tol=1)
    expect_equal(30, vals$J, tol=0.10)
    expect_equal(10, vals$K, tol=0.10)

    expect_equal( warning_pattern(vals), c(TRUE, FALSE,FALSE) )

})



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d3.2_m3rr2rc -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


test_that("testing of d3.2_m3rr2rc", {

    if ( FALSE ) {
        set.seed( 245444 )

        pp1 <- pump_power(
            design = "d3.2_m3rr2rc",
            MTP = 'Holm',
            nbar = 50,
            K = 10,
            J = 30,
            M = 3,
            MDES = rep(0.125, 3),
            Tbar = 0.5, alpha = 0.05,
            numCovar.1 = 1, numCovar.2 = 1,
            R2.1 = 0.1, R2.2 = 0.1,
            ICC.2 = 0.2, ICC.3 = 0.2,
            omega.2 = 0, omega.3 = 0.1, rho = 0.5, tnum = 100000)
        pp1$D1indiv[2]
    }
    pp_power = 0.19405

    vals <- test_sample_triad(target_power = pp_power,
                              nbar = 50, J = 30, K = 10,
                              seed = 30033303,
                              design = "d3.2_m3rr2rc",
                              MTP = 'Holm',
                              power.definition = 'D1indiv',
                              M = 3,
                              MDES = 0.125,
                              Tbar = 0.5, alpha = 0.05,
                              numCovar.1 = 1, numCovar.2 = 1,
                              R2.1 = 0.1, R2.2 = 0.1,
                              ICC.2 = 0.2, ICC.3 = 0.2,
                              omega.2 = 0, omega.3 = 0.1, rho = 0.5 )
    vals[1:3]

    expect_equal(vals$K, 10, tol = 0.1)
    expect_equal(vals$J, 30, tol = 0.1)
    expect_equal(vals$nbar, 50, tol = 1)
    #plot_power_search(vals$nbarrun)

})



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------ d3.3_m3rc2rc -------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# does not converge for nbar

test_that("testing of d3.3_m3rc2rc", {

    if ( FALSE ) {

        set.seed(2344)

        pp1 <- pump_power(
            design = "d3.3_m3rc2rc",
            MTP = 'Holm',
            nbar = 50,
            K = 20,
            J = 40,
            M = 3,
            MDES = rep(0.25, 3),
            Tbar = 0.5, alpha = 0.05,
            numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
            R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
            ICC.2 = 0.1, ICC.3 = 0.1,
            omega.2 = 0, omega.3 = 0, rho = 0.5)
        pp1
        pp1$D1indiv[2]
    }
    pp_power = 0.2594

    vals <- test_sample_triad( target_power = pp_power,
                               nbar = 50, K = 20, J = 40,
                               seed = 4053443,
                               design = "d3.3_m3rc2rc",
                               MTP = 'Holm',
                               power.definition = 'D1indiv',
                               M = 3,
                               MDES = 0.25,
                               Tbar = 0.5, alpha = 0.05,
                               numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
                               R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
                               ICC.2 = 0.1, ICC.3 = 0.1,
                               omega.2 = 0, omega.3 = 0, rho = 0.5)
    vals[1:3]

    expect_equal( 20, vals$K, tol = 0.1)
    expect_equal( 40, vals$J, tol = 0.50)
    expect_true( !is.na( vals$nbar ) )
    #expect_equal( 50, vals$nbar, tol = 0.1)

    #plot_power_search(vals$nbarrun)

    # converges but is very flat
    set.seed( 245444 )
    J1 <- expect_warning(pump_sample(
        design = "d3.3_m3rc2rc",
        typesample = 'J',
        MTP = 'Holm',
        target.power =pp_power,
        power.definition = 'D1indiv',
        K = 20,
        nbar = 50,
        M = 3,
        MDES = 0.25,
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5))
    J1
    expect_true(!is.na(J1$`Sample.size`))
    expect_equal(40, J1$`Sample.size`, tol = 0.1)


    set.seed( 245444 )
    J3 <- expect_warning(pump_sample(
        design = "d3.3_m3rc2rc",
        typesample = 'J',
        MTP = 'Holm',
        target.power =pp_power,
        power.definition = 'D1indiv',
        K = 20,
        nbar = 50,
        M = 3,
        MDES = 0.25,
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5,
        start.tnum = 2000, max.tnum = 4000,
        tol = 0.005, max_sample_size_JK = 80))
    J3
    expect_true(!is.na(J3$`Sample.size`))


    set.seed( 245444 )
    nbar1 <- expect_warning(pump_sample(
        design = "d3.3_m3rc2rc",
        power.definition = 'D1indiv',
        target.power =pp_power,
        typesample = 'nbar',
        MTP = 'Holm',
        K = 20,
        J = 40,
        M = 3,
        MDES = 0.25,
        Tbar = 0.5, alpha = 0.05,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5))
    nbar1
    # plot_power_search(nbar1)

})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------ lower limit -----
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that( "testing of lower limit", {
    # This should hit lower limit (too powerful, want J < 3).
    set.seed( 24553453 )
    expect_warning( pp <- pump_sample(    design = "d3.2_m3ff2rc",
                                          typesample = "J",
                                          MTP = "Holm",
                                          MDES = 0.12,
                                          target.power = 0.50,
                                          power.definition = "min1",
                                          tol = 0.01,
                                          M = 5,
                                          K = 7, # number RA blocks
                                          nbar = 58,
                                          Tbar = 0.50, # prop Tx
                                          alpha = 0.15, # significance level
                                          numCovar.1 = 1, numCovar.2 = 1,
                                          R2.1 = 0.1, R2.2 = 0.7,
                                          ICC.2 = 0.05, ICC.3 = 0.9,
                                          rho = 0.4, # how correlated outcomes are
                                          max.tnum = 200 ) )
    pp
    expect_true( !is.null( pp ) )
    expect_true( pp$`min1 power` > 0.50 )
    expect_true( pp$`Sample.size` == 3 )

} )
