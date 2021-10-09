# library( pum )
# library( testthat )

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
                    tnum = 200,
                    long.table = TRUE
  )
  pp
  expect_true( is.pumpresult( pp ) )

  expect_equal( dim( pp ), c(7,3) )

  expect_true( is.na( pp$None[[5]] ) )
  expect_true( pp$Bonferroni[4] < pp$Bonferroni[5] )
  expect_true( pp$Bonferroni[4] > pp$Bonferroni[7] )
  expect_true( all ( ( pp$None >= pp$Bonferroni )[1:4] ) )
})

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

  expect_equal( dim( pp ), c(2,8) )

  expect_true( pp[2,"min1"] >= pp[2,"D1indiv"] )
})

test_that("having no covariates is fine", {
  pp <- pump_power(   design = "d2.1_m2fc",
                      MTP = "Bonferroni",
                      MDES = rep( 0.10, 3 ),
                      M = 3,
                      J = 3, # number of schools/block
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      ICC.2 = 0.05,
                      rho = 0.4, tnum = 200
  )

  expect_equal( dim( pp ), c(2,8) )
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
  expect_equal( dim( pp ), c(3,8) )

})

test_that("M = 1 runs successfully", {

  pp <- pump_power(   design = "d2.1_m2fc",
                      MTP = "None",
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

test_that("J = 1 runs successfully", {

    expect_warning(pp <- pump_power(   design = "d2.1_m2fc",
                        MTP = "None",
                        MDES = 0.10,
                        M = 1,
                        J = 1, # number of schools/block
                        nbar = 258,
                        Tbar = 0.50, # prop Tx
                        alpha = 0.05, # significance level
                        numCovar.1 = 5, numCovar.2 = 3,
                        R2.1 = 0.1, R2.2 = 0.7,
                        ICC.2 = 0.05,
                        rho = 0.4, tnum = 200
    ))
    expect_true( nrow( pp ) == 1 )
})

test_that("K = 1 runs successfully", {
    expect_warning(pp <- pump_power( design = "d3.2_m3ff2rc",
                      MTP = "Bonferroni",
                      MDES = rep( 0.10, 3 ),
                      M = 3,
                      J = 10, # number of schools/block
                      K = 1, # number RA blocks
                      nbar = 150,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 0, numCovar.2 = 0,
                      R2.1 = 0, R2.2 = 0,
                      ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0.4, tnum = 10000
    ))
    pp
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
             R2.2 = R2.2, numCovar.2 = 1,
             M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5 )


  expect_true( nrow( pp ) == 2 )

})

test_that("Correct MTP parameter validation.", {

    expect_error(pp <- pump_power(   design = "d2.1_m2fc",
                                       MTP = "None",
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
    expect_true( pp.rhomin$min1[2] > pp.rhomed$min1[2] )
    expect_true( pp.rhomed$min1[2] > pp.rhomax$min1[2]  )

    # complete power is the reverse
    expect_true( pp.rhomin$complete[2] <  pp.rhomed$complete[2] )
    expect_true( pp.rhomed$complete[2] < pp.rhomax$complete[2] )

})

test_that("numZero has expected behavior", {

  pp <- pump_power( design = "d2.2_m2rc",
                    MTP = "Bonferroni",
                    J = 10,
                    M = 5,
                    nbar = 100,
                    MDES = rep( 0.2, 2 ),
                    numZero = 3,
                    Tbar = 0.50,
                    rho = 1)



  expect_error(pp <- pump_power( design = "d2.2_m2rc",
                    MTP = "Bonferroni",
                    J = 10,
                    M = 4,
                    nbar = 100,
                    MDES = rep( 0.2, 2 ),
                    numZero = 3,
                    Tbar = 0.50,
                    rho = 1)
  )

})

test_that("do not report invalid power values", {

  pp <- pump_power( design = "d2.2_m2rc",
                    MTP = "Bonferroni",
                    J = 10,
                    M = 3,
                    nbar = 100,
                    MDES = 0.2,
                    Tbar = 0.50,
                    rho = 0.4)

  expect_true(is.na(pp$min1[1]))
  expect_true(is.na(pp$min2[1]))
  expect_true(is.na(pp$complete[1]))
})

