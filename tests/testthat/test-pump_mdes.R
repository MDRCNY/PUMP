
library( pum )
library( testthat )


test_that("pump_mdes runs for Bonferroni", {
  pmdes <- pump_mdes( design="blocked_i1_2c",
                      MTP = "Bonferroni",
                      nbar = 200, J = 50,
                      power.definition = "min1",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      tnum = 300,
                      rho = 0.4 )

  pmdes
  pmdes
  expect_true( pmdes$mdes.results$`Adjusted MDES` > 0 )
  expect_true( abs( pmdes$mdes.results$`D1indiv power` - 0.80) <  0.01 )

  pp = pump_power( design="blocked_i1_2c",
                   MTP = "Bonferroni",
                   MDES = rep( pmdes$mdes.results$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  pp
  expect_true( abs( pp[2,"min1"] - 0.80 ) <= 0.02)
})





test_that("pump_mdes runs for D1indiv, Holm", {
  pmdes <- pump_mdes( design="blocked_i1_2c",
               MTP = "Holm",
               nbar = 200, J = 50,
               power.definition = "D1indiv",
               M = 3,
               target.power = 0.80, tol = 0.01,
               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
               tnum = 300,
               rho = 0.4 )

  pmdes
  expect_true( pmdes$mdes.results$`Adjusted MDES` > 0 )
  expect_true( abs( pmdes$mdes.results$`D1indiv power` - 0.80) <  0.01 )

  pp = pump_power( design="blocked_i1_2c",
                   MTP = "Holm",
                   MDES = rep( pmdes$mdes.results$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  pp
  expect_true( abs( pp[2,1] - 0.80 ) <= 0.02)
})




