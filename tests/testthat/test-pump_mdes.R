
# library( pum )
# library( testthat )


test_that("pump_mdes runs for Bonferroni", {

  pmdesB <- pump_mdes( design = "d2.1_m2fc",
                       MTP = "Bonferroni",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                       R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                       tnum = 300,
                       rho = 0.4 )
  expect_true( pmdesB$mdes.results$`Adjusted MDES` > 0 )
  expect_true( abs(pmdesB$mdes.results$`D2indiv power` - 0.80) <  0.01 )

  pmdesR <- pump_mdes( design = "d2.1_m2fc",
                       MTP = "None",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                       R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                       tnum = 300,
                       rho = 0.4 )
  expect_true( pmdesR$mdes.results$`Adjusted MDES` > 0 )
  expect_true( abs(pmdesR$mdes.results$`D2indiv power` - 0.80) <  0.01 )
  expect_true( pmdesR$mdes.results$`Adjusted MDES` < pmdesB$mdes.results$`Adjusted MDES`)

  pmdesBmin <- pump_mdes(
                      design = "d2.1_m2fc",
                      MTP = "Bonferroni",
                      nbar = 200, J = 50,
                      power.definition = "min1",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      tnum = 300,
                      rho = 0.4 )

  expect_true( pmdesBmin$mdes.results$`Adjusted MDES` < pmdesR$mdes.results$`Adjusted MDES` )
  expect_true( abs( pmdesBmin$mdes.results$`min1 power` - 0.80) <  0.01 )

  pmdesBcomp <- pump_mdes(
                      design = "d2.1_m2fc",
                      MTP = "Bonferroni",
                      nbar = 200, J = 50,
                      power.definition = "complete",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      tnum = 300,
                      rho = 0.4 )
  expect_true( pmdesBcomp$mdes.results$`Adjusted MDES` > pmdesB$mdes.results$`Adjusted MDES` )
  expect_true( abs( pmdesBcomp$mdes.results$`complete power` - 0.80) <  0.01 )

  ppBcomp <- pump_power(
                   design = "d2.1_m2fc",
                   MTP = "Bonferroni",
                   MDES = rep( pmdesBcomp$mdes.results$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )

  expect_true( abs( ppBcomp[2,"complete"] - 0.80 ) <= 0.02)


  ppBmin <- pump_power(
                   design = "d2.1_m2fc",
                   MTP = "Bonferroni",
                   MDES = rep( pmdesBmin$mdes.results$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  expect_true( abs( ppBmin[2,"min1"] - 0.80 ) <= 0.02)

})


test_that("pump_mdes runs for D1indiv, Holm", {
  pmdes <- pump_mdes(
               design = "d2.1_m2fc",
               MTP = "Holm",
               nbar = 200, J = 50,
               power.definition = "D1indiv",
               M = 3,
               target.power = 0.80, tol = 0.01,
               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
               tnum = 300,
               rho = 0.4 )

  pmdes
  expect_true( pmdes$mdes.results$`Adjusted MDES` > 0 )
  expect_true( abs( pmdes$mdes.results$`D1indiv power` - 0.80) <  0.01 )

  pp = pump_power( design = "d2.1_m2fc",
                   MTP = "Holm",
                   MDES = rep( pmdes$mdes.results$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  pp
  expect_true( abs( pp[2,1] - 0.80 ) <= 0.02)
})




