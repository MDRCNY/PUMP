# library( pum )
# library( testthat )

test_that("skipping level three inputs for level 2 works", {

  expect_warning(pp <- pump_power(   design = "d2.1_m2fc",
                                     MTP = "Bonferroni",
                                     MDES = rep( 0.10, 3 ),
                                     M = 3,
                                     J = 3, # number of schools/block
                                     K = 21, # number RA blocks
                                     nbar = 258,
                                     Tbar = 0.50, # prop Tx
                                     alpha = 0.05, # significance level
                                     numCovar.1 = 5, numCovar.2 = 3,
                                     R2.1 = 0.1, R2.2 = 0.7,
                                     ICC.2 = 0.05,
                                     rho = 0.4, tnum = 200
  ))

  expect_equal( dim( pp ), c(2,7) )

  expect_true( all( pp[,"min1"] >= pp[,"D1indiv"] ) )
})


test_that("pump_power works without crashing", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = "Bonferroni",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 21, # number RA blocks
                    nbar = 258,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, # how correlated outcomes are
                    tnum = 200
              )

  expect_equal( dim( pp ), c(2,7) )

  expect_true( all( pp[,"min1"] >= pp[,"D1indiv"] ) )
})



test_that("pump_power works with multiple MTP", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = c( "Holm", "Bonferroni" ),
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 21, # number RA blocks
                    nbar = 258,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.15, # significance level
                    numCovar.1 = 1, numCovar.2 = 1,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.9,
                    rho = 0.4, tnum = 200
  ) # how correlated outcomes are
  expect_equal( dim( pp ), c(3,7) )

})

test_that("M = 1 runs successfully", {

  pp <- pump_power(   design = "d2.1_m2fc",
                      MTP = "none",
                      MDES = 0.10,
                      M = 1,
                      J = 3, # number of schools/block
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5, numCovar.2 = 3,
                      R2.1 = 0.1, R2.2 = 0.7,
                      ICC.2 = 0.05,
                      rho = 0.4, tnum = 200
  )
  expect_true( nrow( pp ) == 1 )
})

test_that("pump_power_grid works", {
  pp <- pump_power_grid(    design = "d3.2_m3ff2rc",
                            MTP = "Holm",
                            MDES = rep( 0.10, 3 ),
                            M = 3,
                            J = c( 3, 5, 9), # number of schools/block
                            K = 21, # number RA blocks
                            nbar = 258,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200
  )
  pp
  expect_equal( nrow(pp), 3 * 2 )

})

test_that("unblocked designs", {

  pp <- pump_power(   design = "d1.1_m2cc",
                      MTP = "Bonferroni",
                      MDES = rep( 0.50, 3 ),
                      M = 3,
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5, numCovar.2 = 3,
                      R2.1 = 0.1, R2.2 = 0.7,
                      ICC.2 = 0.05,
                      rho = 0.4, tnum = 200
  )


  ES = log( 2 ) / 0.66
  ES
  R2.2 = 0.6102
  pump_power(design = "d1.1_m2cc", MTP = "Holm", MDES = ES,
             R2.2 = R2.2,
             M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5 )


  expect_true( nrow( pp ) == 2 )

})

test_that("Correct MTP parameter validation. ", {

    expect_error(pp <- pump_power(   design = "d2.1_m2fc",
                                       MTP = "rawp",
                                       MDES = rep( 0.10, 3 ),
                                       M = 3,
                                       J = 3, # number of schools/block
                                       K = 21, # number RA blocks
                                       nbar = 258,
                                       Tbar = 0.50, # prop Tx
                                       alpha = 0.05, # significance level
                                       numCovar.1 = 5, numCovar.2 = 3,
                                       R2.1 = 0.1, R2.2 = 0.7,
                                       ICC.2 = 0.05,
                                       rho = 0.4, tnum = 200
    ))

    # no MTP provided
    expect_error(pp <- pump_power( design = "d2.1_m2fc",
                                     MDES = 0.1,
                                     M = 3,
                                     J = 3, # number of schools/block
                                     nbar = 258,
                                     Tbar = 0.50, # prop Tx
                                     alpha = 0.05, # significance level
                                     numCovar.1 = 5, numCovar.2 = 3,
                                     R2.1 = 0.1, R2.2 = 0.7,
                                     ICC.2 = 0.05,
                                     rho = 0.4, tnum = 200
    ))

    # no MTP provided with single outcome is fine.
    pp <- pump_power( design = "d2.1_m2fc",
                                   MDES = 0.1,
                                   M = 1,
                                   J = 3, # number of schools/block
                                   nbar = 258,
                                   Tbar = 0.50, # prop Tx
                                   alpha = 0.05, # significance level
                                   numCovar.1 = 5, numCovar.2 = 3,
                                   R2.1 = 0.1, R2.2 = 0.7,
                                   ICC.2 = 0.05,
                                   rho = 0.4, tnum = 200
    )
    expect_true( nrow( pp ) == 1 )
})


test_that("different correlations", {

    pp.rhomin <- pump_power( design = "d2.2_m2rc",
                             MTP = "Bonferroni",
                             J = 10,
                             M = 4,
                             nbar = 100,
                             MDES = rep( 0.2, 4 ),
                             Tbar = 0.50, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.5, ICC.2 = 0.05,
                             rho = 0)

    pp.rhomed <- pump_power(   design = "d2.2_m2rc",
                               MTP = "Bonferroni",
                               J = 10,
                               M = 4,
                               nbar = 100,
                               MDES = rep( 0.2, 4 ),
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.5, ICC.2 = 0.05,
                               rho = 0.4)

    pp.rhomax <- pump_power(   design = "d2.2_m2rc",
                               MTP = "Bonferroni",
                               J = 10,
                               M = 4,
                               nbar = 100,
                               MDES = rep( 0.2, 4 ),
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.5, ICC.2 = 0.05,
                               rho = 1)

    pp.rhomin
    pp.rhomed
    pp.rhomax
    # lower correlation means higher min1 power (more chances to hit one out of the park)
    expect_true( all( pp.rhomin$min1 > pp.rhomed$min1 ) )
    expect_true( all( pp.rhomed$min1 > pp.rhomax$min1 ) )

    # complete power is the reverse
    expect_true( all( pp.rhomin$complete <  pp.rhomed$complete ) )
    expect_true( all( pp.rhomed$complete < pp.rhomax$complete ) )

})

