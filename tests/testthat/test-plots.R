# library(testthat)

test_that("Single Scenario plot works", {

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
  ss.plot <- plot(pp)

  expect_true(!is.null(ss.plot))

})

test_that("Grid plot works", {

  grid <- pump_power_grid( design = "d3.2_m3ff2rc",
                           MTP = "Bonferroni",
                           MDES = 0.1,
                           M = 3,
                           J = 3, # number of schools/block
                           K = 21, # number RA blocks
                           nbar = 258,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, numCovar.2 = 3,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = 0.3,
                           ICC.3 = seq( 0, 0.45, 0.15 ),
                           rho = 0.4,
                           tnum = 100,
                           long.table = TRUE)
  grid.plot <- plot(grid, power.definition = 'min1', var.vary = 'ICC.3')

  expect_true(!is.null(grid.plot))


  grid <- pump_power_grid( design = "d3.2_m3ff2rc",
                           MTP = "Bonferroni",
                           MDES = 0.1,
                           M = 3,
                           J = 3, # number of schools/block
                           K = 21, # number RA blocks
                           nbar = 258,
                           Tbar = 0.50, # prop Tx
                           alpha = 0.05, # significance level
                           numCovar.1 = 5, numCovar.2 = 3,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = seq( 0, 0.45, 0.15 ),
                           ICC.3 = seq( 0, 0.45, 0.15 ),
                           rho = 0.4,
                           tnum = 100,
                           long.table = TRUE)

  expect_error(grid.plot <- plot(grid, power.definition = 'min1', var.vary = 'ICC.3') )
})

