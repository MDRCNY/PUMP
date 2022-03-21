# library( PUMP )
# library( testthat )

default.tnum <- 1000

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ----- three level models ------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# --------    d3.1_m3rr2rr    --------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

skip_on_cran()

test_that("testing of d3.1_m3rr2rr one-tailed", {

    if ( FALSE ) {

        set.seed( 524235326 )
        pp1 <- pump_power(
            d_m = "d3.1_m3rr2rr",
            MTP = 'HO',
            nbar = 50,
            K = 15,
            J = 30,
            M = 3,
            MDES = rep(0.125, 3),
            Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
            numCovar.1 = 1, numCovar.2 = 1,
            R2.1 = 0.1, R2.2 = 0.1,
            ICC.2 = 0.2, ICC.3 = 0.2,
            omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
            tnum = 100000)
        pp1
        pp_power <- pp1$D1indiv[2]
    }

    pp_power <- 0.88166
    
    vals <- test_sample_triad( 
      target_power = pp_power, nbar = 50, J = 30, K = 15,
      seed = 4224423,
      d_m = "d3.1_m3rr2rr",
      MTP = 'HO',
      power.definition = 'D1indiv',
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.1,
      ICC.2 = 0.2, ICC.3 = 0.2,
      omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
      tnum = default.tnum )
    
    
    vals
    
    expect_equal(15, vals$K, tolerance = 0.1)
    # expect_equal(30, vals$J, tolerance = 0.1)
    expect_equal( warning_pattern(vals), c(TRUE, TRUE, FALSE) )
    
    
    
    set.seed( 44040422 )
    mdes1 <-  pump_mdes(
      d_m = "d3.1_m3rr2rr",
      MTP = 'HO',
      power.definition = 'D1indiv',
      target.power = pp_power,
      J = 30,
      K = 15,
      nbar = 50,
      M = 3,
      Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.1,
      ICC.2 = 0.2, ICC.3 = 0.2,
      omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
      tnum = default.tnum )
    
    expect_equal(0.125, mdes1$Adjusted.MDES, tolerance = 0.1)

    # if we go below the true value, we get the wrong number since it is so flat
    set.seed( 524235325 )

    nbar2 <- pump_sample(
        d_m = "d3.1_m3rr2rr",
        typesample = 'nbar',
        MTP = 'HO',
        target.power = 0.66682,
        power.definition = 'D1indiv',
        K = 15,
        J = 30,
        M = 3,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
        numCovar.1 = 1, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.1,
        ICC.2 = 0.2, ICC.3 = 0.2,
        omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
        tnum = default.tnum,
        max_sample_size_nbar = 40 )
    expect_true(nbar2$`Sample.size` < 40 )
})



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------ d3.2_m3ff2rc ------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


test_that("testing of d3.2_m3ff2rc two-tailed", {
    
    if ( FALSE ) {

      set.seed( 245444 )
      pp1 <- pump_power(
          d_m = "d3.2_m3ff2rc",
          MTP = 'HO',
          nbar = 50,
          J = 30,
          K = 10,
          M = 5,
          MDES = 0.125,
          Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
          numCovar.1 = 1, numCovar.2 = 1,
          R2.1 = 0.1, R2.2 = 0.1,
          ICC.2 = 0.2, ICC.3 = 0.2,
          omega.2 = 0, omega.3 = 0.1, rho = 0.5, tnum = 100000)
      pp1
      pp_power <- pp1$min2[2]
    }
    pp_power <- 0.64854

    set.seed( 245444 )
    vals <- test_sample_triad( pp_power, nbar = 50, J = 30, K = 10,
                               seed = 4224422,
                               d_m = "d3.2_m3ff2rc",
                               MTP = 'HO',
                               power.definition = 'min2',
                               M = 5,
                               MDES = 0.125,
                               Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
                               numCovar.1 = 1, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.1,
                               ICC.2 = 0.2, ICC.3 = 0.2,
                               omega.2 = 0, omega.3 = 0.1, rho = 0.5,
                               tnum = default.tnum )
    vals[1:3]

    # nbar ends up not converging
    #expect_equal(50, vals$nbar, tol=0.5 )
    expect_equal(30, vals$J, tol = 0.10)
    expect_equal(10, vals$K, tol = 0.10)

    expect_equal( warning_pattern(vals), c(TRUE, FALSE, FALSE) )

    # nbar converges but is flat
    set.seed( 245444 )
    pp_power
    nbar1 <- expect_warning( pump_sample(
      d_m = "d3.2_m3ff2rc",
      MTP = 'HO',
      power.definition = 'min2',
      typesample = 'nbar',
      target.power = pp_power,
      M = 5,
      J = 30, K = 10,
      MDES = 0.125,
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.1,
      ICC.2 = 0.2, ICC.3 = 0.2,
      omega.2 = 0, omega.3 = 0.1, rho = 0.5, max.steps = 40,
      max_sample_size_nbar = 1000, tnum = 4000))
    nbar1
    # plot_power_search(nbar1)
    expect_equal(50, nbar1$`Sample.size`, tolerance = 0.4)
    expect_true( attr( nbar1, "flat" ) )
})



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d3.2_m3rr2rc -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


test_that("testing of d3.2_m3rr2rc one tailed", {
    
    if ( FALSE ) {
        set.seed( 245444 )

        pp1 <- pump_power(
            d_m = "d3.2_m3rr2rc",
            MTP = 'HO',
            nbar = 50,
            K = 10,
            J = 30,
            M = 3,
            MDES = rep(0.125, 3),
            Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
            numCovar.1 = 1, numCovar.2 = 1,
            R2.1 = 0.1, R2.2 = 0.1,
            ICC.2 = 0.2, ICC.3 = 0.2,
            omega.2 = 0, omega.3 = 0.1, rho = 0.5,
            tnum = 100000)
        pp_power <- pp1$D1indiv[2]
    }
    pp_power <- 0.33201

    vals <- test_sample_triad(target_power = pp_power,
                              nbar = 50, J = 30, K = 10,
                              seed = 30033303,
                              d_m = "d3.2_m3rr2rc",
                              MTP = 'HO',
                              power.definition = 'D1indiv',
                              M = 3,
                              MDES = 0.125,
                              Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
                              numCovar.1 = 1, numCovar.2 = 1,
                              R2.1 = 0.1, R2.2 = 0.1,
                              ICC.2 = 0.2, ICC.3 = 0.2,
                              omega.2 = 0, omega.3 = 0.1, rho = 0.5,
                              tnum = default.tnum )
    vals[1:3]

    # nbar is flat!
    expect_equal(vals$K, 10, tol = 0.1)
    expect_equal(vals$J, 30, tol = 0.1)
    #expect_equal(vals$nbar, 50, tol = 0.4)

})



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------ d3.3_m3rc2rc -------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d3.3_m3rc2rc two tailed", {
    
    set.seed(2344)

    if ( FALSE ) {

        set.seed(2344)

        pp1 <- pump_power(
            d_m = "d3.3_m3rc2rc",
            MTP = 'HO',
            nbar = 50,
            K = 20,
            J = 40,
            M = 3,
            MDES = rep(0.25, 3),
            Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
            numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
            R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
            ICC.2 = 0.1, ICC.3 = 0.1,
            omega.2 = 0, omega.3 = 0, rho = 0.5,
            tnum = 100000)
        pp1
        pp_power <- pp1$D1indiv[2]
    }

    pp_power <- 0.25873

    vals <- test_sample_triad( target_power = pp_power,
                               nbar = 50, K = 20, J = 40,
                               seed = 4053443,
                               d_m = "d3.3_m3rc2rc",
                               MTP = 'HO',
                               power.definition = 'D1indiv',
                               M = 3,
                               MDES = 0.25,
                               Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
                               numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
                               R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
                               ICC.2 = 0.1, ICC.3 = 0.1,
                               omega.2 = 0, omega.3 = 0, rho = 0.5,
                               tnum = default.tnum )
    vals[1:3]

    expect_equal( 20, vals$K, tol = 0.1)
    #expect_equal( 40, vals$J, tol = 0.50)
    #expect_true( is.na( vals$nbar ) )
    expect_true( length( vals$Kwarn ) == 0 )
    expect_true( length( vals$Jwarn ) > 0 )
    expect_true( length( vals$nbarwarn ) > 0 )
    

    # converges but is relatively flat
    set.seed( 245444 )
    J1 <- expect_warning(pump_sample(
        d_m = "d3.3_m3rc2rc",
        typesample = 'J',
        MTP = 'HO',
        target.power = pp_power,
        power.definition = 'D1indiv',
        K = 20,
        nbar = 50,
        M = 3,
        MDES = 0.25,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5,
        tnum = default.tnum ))
    J1
    expect_true(!is.na(J1$`Sample.size`))
    expect_true( attr(J1, "flat") ) 
    
    # very flat!
    set.seed( 245444 )
    J2 <- expect_warning(pump_sample(
        d_m = "d3.3_m3rc2rc",
        typesample = 'J',
        MTP = 'HO',
        target.power = pp_power,
        power.definition = 'D1indiv',
        K = 20,
        nbar = 50,
        M = 3,
        MDES = 0.25,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
        R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
        ICC.2 = 0.1, ICC.3 = 0.1,
        omega.2 = 0, omega.3 = 0, rho = 0.5,
        tnum = 10000, final.tnum = 10000,
        tol = 0.005, max_sample_size_JK = 80))
    J2
    expect_true(!is.na(J2$`Sample.size`))
    expect_true( attr(J2, "flat") ) 
    
})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------ lower limit -----
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that( "testing of lower limit", {
    
    # This should hit lower limit (too powerful, want J < 3).
    set.seed( 24553453 )
    expect_warning( pp <- pump_sample(    d_m = "d3.2_m3ff2rc",
                                          typesample = "J",
                                          MTP = "HO",
                                          MDES = 0.12,
                                          target.power = 0.50,
                                          power.definition = "min1",
                                          tol = 0.02,
                                          M = 5,
                                          K = 7, # number RA blocks
                                          nbar = 58,
                                          Tbar = 0.50, # prop Tx
                                          alpha = 0.15, two.tailed = TRUE, # significance level
                                          numCovar.1 = 1, numCovar.2 = 1,
                                          R2.1 = 0.1, R2.2 = 0.7,
                                          ICC.2 = 0.05, ICC.3 = 0.9,
                                          rho = 0.4, # how correlated outcomes are
                                          tnum = 200 ) )
    pp
    expect_true( !is.null( pp ) )
    expect_true( pp$`min1 power` > 0.50 )
    expect_true( pp$`Sample.size` == 3 )

} )

