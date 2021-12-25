# library( PUMP )
# library( testthat )

small.B <- 1000
default.tol <- 0.05

test_that("WY-SS results are stable for design: d2.1_m2fc (constant)", {

  pp1 <- pump_power( design = "d2.1_m2fc",
                       MTP = "WY-SS",
                       MDES = rep( 0.10, 3 ),
                       M = 3,
                       J = 10, # number of schools/block
                       nbar = 300,
                       Tbar = 0.50, # prop Tx
                       alpha = 0.05, # significance level
                       numCovar.1 = 5,
                       R2.1 = 0.1, 
                       ICC.2 = 0.05,
                       rho = 0.4, tnum = 1000,
                       B = small.B
  )

  pp2 <- pump_power( design = "d2.1_m2fc",
                          MTP = "WY-SS",
                          MDES = rep( 0.10, 3 ),
                          M = 3,
                          J = 10, # number of schools/block
                          nbar = 300,
                          Tbar = 0.50, # prop Tx
                          alpha = 0.05, # significance level
                          numCovar.1 = 5,
                          R2.1 = 0.1,
                          ICC.2 = 0.05,
                          rho = 0.4, tnum = 1000,
                          B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})

test_that("WY-SD results are stable for design: d2.1_m2fc (constant)", {

  pp1 <- pump_power( design = "d2.1_m2fc",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 10, # number of schools/block
                     nbar = 300,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5,
                     R2.1 = 0.1,
                     ICC.2 = 0.05,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp2 <- pump_power( design = "d2.1_m2fc",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 10, # number of schools/block
                     nbar = 300,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5,
                     R2.1 = 0.1,
                     ICC.2 = 0.05,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})

test_that("WY-SS results are stable for design: d3.1_m3rr2rr (random)", {

  pp1 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SS",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 10, # number of schools/block
                     K = 10,
                     nbar = 300,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.1, R2.2 = 0.4, R2.3 = 0.3,
                     ICC.2 = 0.05, ICC.3 = 0.2,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )


  pp2 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SS",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 10, # number of schools/block
                     K = 10,
                     nbar = 300,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.1, R2.2 = 0.4, R2.3 = 0.3,
                     ICC.2 = 0.05, ICC.3 = 0.2,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})

test_that("WY-SD results are stable for design: d3.1_m3rr2rr (random)", {

  pp1 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 10, # number of schools/block
                     K = 10,
                     nbar = 300,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.1, R2.2 = 0.4, R2.3 = 0.3,
                     ICC.2 = 0.05, ICC.3 = 0.2,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )


  pp2 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 10, # number of schools/block
                     K = 10,
                     nbar = 300,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.1, R2.2 = 0.4, R2.3 = 0.3,
                     ICC.2 = 0.05, ICC.3 = 0.2,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})


test_that("WY-SS results are stable for design: d3.1_m3rr2rr (random), small number of blocks/clusters", {

  pp1 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SS",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 5, # number of schools/block
                     K = 5,
                     nbar = 500,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.5, R2.2 = 0.5, R2.3 = 0.5,
                     ICC.2 = 0.05, ICC.3 = 0.05,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )


  pp2 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SS",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 5, # number of schools/block
                     K = 5,
                     nbar = 500,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.5, R2.2 = 0.5, R2.3 = 0.5,
                     ICC.2 = 0.05, ICC.3 = 0.05,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})


test_that("WY-SD results are stable for design: d3.1_m3rr2rr (random), small number of blocks/clusters", {

  pp1 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 5, # number of schools/block
                     K = 5,
                     nbar = 500,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.5, R2.2 = 0.5, R2.3 = 0.5,
                     ICC.2 = 0.05, ICC.3 = 0.05,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )


  pp2 <- pump_power( design = "d3.1_m3rr2rr",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 5, # number of schools/block
                     K = 5,
                     nbar = 500,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 5, numCovar.2 = 3, numCovar.3 = 1,
                     R2.1 = 0.5, R2.2 = 0.5, R2.3 = 0.5,
                     ICC.2 = 0.05, ICC.3 = 0.05,
                     omega.3 = 0.1, omega.2 = 0.1,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})


test_that("WY-SS results are stable for design: d2.2_m2rc", {

  pp1 <- pump_power( design = "d2.2_m2rc",
                     MTP = "WY-SS",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 20, # number of schools/block
                     nbar = 1000,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 1, numCovar.2 = 1,
                     R2.1 = 0.1, R2.2 = 0.1,
                     ICC.2 = 0.01,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )


  pp2 <- pump_power( design = "d2.2_m2rc",
                     MTP = "WY-SS",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 20, # number of schools/block
                     nbar = 1000,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 1, numCovar.2 = 1,
                     R2.1 = 0.1, R2.2 = 0.1,
                     ICC.2 = 0.01,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})


test_that("WY-SD results are stable for design: d2.2_m2rc", {

  pp1 <- pump_power( design = "d2.2_m2rc",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 20, # number of schools/block
                     nbar = 1000,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 1, numCovar.2 = 1,
                     R2.1 = 0.1, R2.2 = 0.1,
                     ICC.2 = 0.01,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )


  pp2 <- pump_power( design = "d2.2_m2rc",
                     MTP = "WY-SD",
                     MDES = rep( 0.10, 3 ),
                     M = 3,
                     J = 20, # number of schools/block
                     nbar = 1000,
                     Tbar = 0.50, # prop Tx
                     alpha = 0.05, # significance level
                     numCovar.1 = 1, numCovar.2 = 1,
                     R2.1 = 0.1, R2.2 = 0.1,
                     ICC.2 = 0.01,
                     rho = 0.4, tnum = 1000,
                     B = small.B
  )

  pp1
  pp2

  expect_equal( pp1[2,2], pp2[2,2], tol = default.tol)
})
