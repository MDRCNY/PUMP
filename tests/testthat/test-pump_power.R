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
                    numCovar.1 = 0, numCovar.2 = 0,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.9,
                    rho = 0.4, tnum = 200
  ) # how correlated outcomes are
  expect_equal( dim( pp ), c(3,7) )

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
                            numCovar.1 = 0, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200
  )
  pp
  expect_equal( nrow(pp), 3 * 2 )

})


test_that("parameters that result in 80% power for raw D1indiv", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = "Bonferroni",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 22, # number RA blocks
                    nbar = 150,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, tnum = 10000
  )
  pp
  expect_equal( pp[1,1], 0.8, tol = 0.05)
})


test_that("parameters that result in 80% power for Bonferroni D1indiv", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = "Bonferroni",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 4, # number of schools/block
                    K = 17, # number RA blocks
                    nbar = 160,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.3, R2.2 = 0.75,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, tnum = 10000
  )

  expect_equal( pp[2,1], 0.8, tol = 0.05)
})


test_that("parameters that result in 80% power for Bonferroni min1", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = "Bonferroni",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 17, # number RA blocks
                    nbar = 100,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.25, R2.2 = 0.75,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, tnum = 10000
  )

  expect_equal( pp[2,5], 0.8, tol = 0.05)
})

test_that("parameters that result in 80% power for Bonferroni min2", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = "Bonferroni",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 4, # number of schools/block
                    K = 16, # number RA blocks
                    nbar = 150,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.25, R2.2 = 0.75,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, tnum = 10000
  )

  expect_equal( pp[2,6], 0.8, tol = 0.05)
})

test_that("parameters that result in 80% power for Bonferroni complete", {
  pp <- pump_power( design = "d3.2_m3ff2rc",
                    MTP = "Bonferroni",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 5, # number of schools/block
                    K = 18, # number RA blocks
                    nbar = 170,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.25, R2.2 = 0.75,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, tnum = 10000
  )
  expect_equal( pp[2,7], 0.8, tol = 0.05)
})

