
# library( pum )
# library( testthat )


test_that("pump_mdes runs for Bonferroni", {
  set.seed( 2424424 )

  pmdesB <- pump_mdes( design = "d2.1_m2fc",
                       MTP = "Bonferroni",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                       R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                       max.tnum = 300,
                       rho = 0.4, just.result.table = FALSE )
  pmdesB
  expect_true( pmdesB$`Adjusted MDES` > 0 )
  expect_true( abs(pmdesB$`D2indiv power` - 0.80) <  0.01 )

  pmdesR <- pump_mdes( design = "d2.1_m2fc",
                       MTP = "None",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                       R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                       max.tnum = 300,
                       rho = 0.4 )
  expect_true( pmdesR$`Adjusted MDES` > 0 )
  expect_true( abs(pmdesR$`D2indiv power` - 0.80) <  0.01 )
  expect_true( pmdesR$`Adjusted MDES` < pmdesB$`Adjusted MDES`)

  set.seed( 14444444 )
  pmdesBmin <- pump_mdes(
                      design = "d2.1_m2fc",
                      MTP = "Bonferroni",
                      nbar = 200, J = 50,
                      power.definition = "min1",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      max.tnum = 1000,
                      rho = 0.4, just.result.table = FALSE )
  pmdesBmin
  expect_true( pmdesBmin$`Adjusted MDES` < pmdesR$`Adjusted MDES` )
  expect_true( abs( pmdesBmin$`min1 power` - 0.80) <  0.01 )

  set.seed( 444224 )
  pmdes_comp <- pump_mdes( design = "d2.1_m2fc",
                           MTP = "Bonferroni",
                           nbar = 200, J = 50,
                           power.definition = "complete",
                           M = 3,
                           target.power = 0.80, tol = 0.01,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                           max.tnum = 300,
                           rho = 0.4, just.result.table = FALSE )
  pmdes_comp
  expect_true( pmdes_comp$`Adjusted MDES` > pmdesB$`Adjusted MDES` )
  expect_true( abs( pmdes_comp$`complete power` - 0.80) <  0.01 )

  ppBcomp <- pump_power(
                   design = "d2.1_m2fc",
                   MTP = "Bonferroni",
                   MDES = rep( pmdes_comp$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )

  ppBcomp
  expect_true( abs( ppBcomp[2,"complete"] - 0.80 ) <= 0.02)

  pmdesBmin$mdes.results
  ppBmin <- pump_power(
                   design = "d2.1_m2fc",
                   MTP = "Bonferroni",
                   MDES = rep( pmdesBmin$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  ppBmin
  expect_true( abs( ppBmin[2,"min1"] - 0.80 ) <= 0.02)

})


test_that("pump_mdes runs for D1indiv, Holm", {

  set.seed( 1010101 )
  pmdes <- pump_mdes( design = "d2.1_m2fc",
                      MTP = "Holm",
                      nbar = 200, J = 50,
                      power.definition = "D1indiv",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      max.tnum = 300,
                      rho = 0.4, just.result.table = FALSE )

  pmdes
  expect_true( pmdes$`Adjusted MDES` > 0 )
  expect_true( abs( pmdes$`D1indiv power` - 0.80) <  0.01 )

  pp = pump_power( design = "d2.1_m2fc",
                   MTP = "Holm",
                   MDES = rep( pmdes$`Adjusted MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  pp
  expect_true( abs( pp[2,1] - 0.80 ) <= 0.02)
})



test_that("pump_mdes runs for d1.1_m2cc", {

  set.seed( 10130103 )
  R2.1 = 0.61
  pmdes <- pump_mdes(design = "d1.1_m2cc", MTP = "Holm",
                     target.power = 0.80, power.definition = "min1", tol = 0.02,
                     R2.1 = R2.1, numCovar.1 = 1, J = 1,
                     max.tnum = 1000,
                     M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5,
                     just.result.table = FALSE )
  pmdes

  ES = pmdes$`Adjusted MDES`
  ppow <- pump_power(design = "d1.1_m2cc", MTP = "Holm", MDES = ES,
             R2.1 = R2.1, numCovar.1 = 1,
             M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5 )
  ppow

  expect_true( !is.null( pmdes ) )
  expect_true( abs( ppow$min1[[2]]- 0.80 ) <= 0.02 )
} )


test_that("No adjustment", {

    pmdes <- pump_mdes( design = "d2.1_m2fc",
                        MTP = "None",
                        nbar = 200, J = 50,
                        power.definition = "D1indiv",
                        M = 3,
                        target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                        max.tnum = 300,
                        rho = 0.4, just.result.table = FALSE )

    expect_error(pmdes <- pump_mdes( design = "d2.1_m2fc",
                        MTP = "None",
                        nbar = 200, J = 50,
                        power.definition = "min2",
                        M = 3,
                        target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                        max.tnum = 300,
                        rho = 0.4, just.result.table = FALSE ))


    pmdes <- pump_mdes( design = "d2.1_m2fc",
                        MTP = "Holm",
                        nbar = 200, J = 50,
                        power.definition = "D1indiv",
                        M = 3,
                        target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                        max.tnum = 300,
                        rho = 0.4, just.result.table = FALSE )
})


