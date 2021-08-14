
test_that("pump_power_grid works", {

  pp <- pump_power_grid(    design = "d3.2_m3ff2rc",
                            MTP = "Holm",
                            MDES = c( 0.05, 0.2 ),
                            M = 5,
                            J = c( 3, 5, 9), # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200,
                            verbose = FALSE
  )
  pp
  expect_equal( nrow(pp), 3 * 2 * 2 )



  pp <- pump_power_grid(    design = "d3.2_m3ff2rc",
                            MTP = "Holm",
                            MDES = c( 0.05, 0.2 ),
                            M = 5,
                            J = 4, # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200,
                            verbose = FALSE
  )
  pp
  expect_equal( nrow(pp), 2 * 2 )



  pp <- pump_power_grid(    design = "d3.2_m3ff2rc",
                            MTP = "Holm",
                            MDES = c( 0.05, 0.2 ),
                            numZero = c(1,2,3),
                            M = 5,
                            J = 4, # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            tnum = 200,
                            verbose = FALSE
  )
  pp
  expect_equal( nrow(pp), 2 * 3 * 2)

})






test_that("pump_mdes_grid works", {

  pp <- pump_mdes_grid(    design = "d3.2_m3ff2rc",
                            MTP = "Holm",
                           target.power = c( 0.50, 0.80 ),
                           power.definition = "D1indiv",
                           tol = 0.05,
                            M = 5,
                            J = c( 3, 9), # number of schools/block
                            K = 7, # number RA blocks
                            nbar = 58,
                            Tbar = 0.50, # prop Tx
                            alpha = 0.15, # significance level
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7,
                            ICC.2 = 0.05, ICC.3 = 0.9,
                            rho = 0.4, # how correlated outcomes are
                            verbose = FALSE, max.tnum = 1000,
  )
  pp
  expect_equal( nrow(pp), 2 * 2)

})







test_that("pump_sample_grid works", {

  pp <- pump_sample_grid(    design = "d3.2_m3ff2rc",
                             typesample = "J",
                           MTP = "Holm",
                           MDES = 0.10,
                           target.power = c( 0.50, 0.80 ),
                           power.definition = "min1",
                           tol = 0.03,
                           M = 5,
                           K = 7, # number RA blocks
                           nbar = 58,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.15, # significance level
                           numCovar.1 = 1, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = 0.25, ICC.3 = 0.25,
                           rho = 0.4, # how correlated outcomes are
                           tnum = 200, verbose = FALSE, max.tnum = 400

  )
  pp
  expect_equal( nrow(pp), 2 )

})


