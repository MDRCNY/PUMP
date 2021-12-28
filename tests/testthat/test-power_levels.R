# library( PUMP )
# library( testthat )

test_that("parameters that result in 80% power for raw D1indiv", {
  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
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
  expect_equal( pp$D1indiv[1], 0.8, tol = 0.05)
})


test_that("parameters that result in 80% power for BF D1indiv", {
  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
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

  expect_equal(pp$D1indiv[2], 0.8, tol = 0.05)
})


test_that("parameters that result in 80% power for BF min1", {
  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
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
  pp
  expect_equal( pp$min1[2], 0.8, tol = 0.05)
})

test_that("parameters that result in 80% power for BF min2", {
  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
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

  expect_equal( pp$min2[2], 0.8, tol = 0.05)
})

test_that("parameters that result in 80% power for BF complete", {
  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
                    MDES = rep( 0.09, 3 ),
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
  expect_equal( pp$complete[2], 0.8, tol = 0.05)
})
