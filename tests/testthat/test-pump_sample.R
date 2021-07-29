
# library( pum )
# library( testthat )


test_that("calc.nbar works", {

  nbar <- pum:::calc.nbar(  design = "d2.2_m2rc",
                            MT = 2.8,
                            MDES = 0.20,
                            J = 5,
                            Tbar = 0.25,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  nbar
  expect_true( is.na( nbar ) )

  nbar <- pum:::calc.nbar(  design = "d2.2_m2rc",
                            MT = 2.8,
                            MDES = 0.20,
                            J = 305,
                            Tbar = 0.25,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  nbar
  expect_true( !is.na( nbar ) )

} )




test_that("calc.J works", {

  J <- pum:::calc.J(  design = "d2.2_m2rc",
                      MT = 2.8,
                      MDES = 0.20,
                      nbar = 200,
                      Tbar = 0.25,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  J
  expect_true( J < 100 )
} )



test_that("pump_sample_raw works", {

  expect_error( calcnbar <- pump_sample_raw(
                            design = "d2.2_m2rc",
                            typesample = "nbar",
                            J = 5,
                            MDES = 0.05, target.power = 0.80, tol = 0.01,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  )

  calcnbar <- pump_sample_raw( design = "d2.2_m2rc",
                               typesample = "nbar",
                               J = 10,
                               MDES = 0.05, target.power = 0.80,
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcnbar
  expect_true( is.na( calcnbar ) )


  calcJ <- pump_sample_raw( design = "d2.1_m2fc",
                            typesample = "J",
                            nbar = 258,
                            MDES = 0.05, target.power = 0.80,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcJ
  expect_true( !is.na( calcJ ) )

  calcn <- pump_sample_raw( design = "d2.1_m2fc",
                            typesample = "nbar",
                            J = calcJ,
                            MDES = 0.05, target.power = 0.80,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcn
  expect_true( abs( calcn - 258 ) < 2 )

  calcJ2 <- pump_sample_raw( design = "d2.1_m2fc",
                             typesample = "J",
                             nbar = calcn,
                             MDES = 0.05, target.power = 0.80,
                             Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                             R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcJ2
  expect_true( abs( calcJ - calcJ2 ) < 1 )

  calcn2 <- pump_sample_raw( design = "d2.1_m2fc",
                             typesample = "nbar",
                             J = calcJ2,
                             MDES = 0.05, target.power = 0.80,
                             Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                             R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcn2
  expect_true( abs( calcn - calcn2 ) < 2 )
})


test_that( "optimize_power solves", {
  op_pow <- pum:::optimize_power(
    MTP = "Holm", nbar=200,
    power.definition="D1indiv",
    design = "d2.1_m2fc", search.type = "J",
    start.low = 56, start.high = 75,
    start.tnum = 200,
    M = 3,
    MDES = 0.05, target.power = 0.80, tol = 0.01,
    Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
    rho = 0.4, max.cum.tnum = 1000, final.tnum = 2000
  )
  op_pow
  expect_true( ncol( op_pow ) == 6 )
  expect_true( all( op_pow$w <= 2000 ) )
  expect_true( max( op_pow$w ) == 2000 )

})


test_that("pump_sample 2 level/2 level", {
  pwr <- pump_sample(   design = "d2.1_m2fc",
                        MTP = "Holm",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                        rho = 0.4 )
  pwr

  p2 = pump_power( design = "d2.1_m2fc",
                   MTP = "Holm",
                   J = pwr$ss.results$`Sample size`,
                   nbar = 200,
                   M = 3,
                   MDES = rep(0.05, 3),
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                   rho = 0.4 )
  p2

  expect_true( abs( p2[ 2, "indiv.mean" ] - 0.80) <= 0.01 )
} )


test_that("sample search when one end is missing", {

  pow_ref <- pump_power( design = "d2.2_m2rc",
                         MTP = "Holm",
                         M = 4,
                         J = 10,
                         nbar = 10000,
                         MDES = rep( 0.40, 4 ),
                         Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                         R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                         rho = 0.2)

  pow_ref

  expect_warning( calcnbar <- pump_sample(
                    design = "d2.2_m2rc",
                    typesample = "nbar",
                    power.definition = "min1",
                    MTP = "Holm",
                    M = 4,
                    J = 10,
                    MDES = 0.40, target.power = 0.80, tol = 0.01,
                    Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                    rho = 0.2) )

  calcnbar
  expect_true( !is.na( calcnbar$ss.results$`Sample size` ) )

  # Now an infeasible calculation where the correlation makes min1 not able to
  # achieve power, even though independence would.
  calcnbar <- pump_sample( design = "d2.2_m2rc",
                           typesample = "nbar",
                           power.definition = "min1",
                           MTP = "Holm",
                           M = 4,
                           J = 10,
                           MDES = 0.39, target.power = 0.80, tol = 0.01,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                           rho = 0.2)
  calcnbar
  expect_true( is.na( calcnbar$ss.results$`Sample size` ) )
})






test_that("further testing of d2.2_m2rc", {


  expect_warning( traw <- pump_sample_raw(
               design = "d2.2_m2rc",
               typesample = "J",
               nbar = 1000,
               MDES = 0.40, target.power = 0.80,
               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 ) )

  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 10,
                           MDES = 0.40, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 10,
                           MDES = 0.01, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 1000,
                           MDES = 0.001, target.power = 0.99,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw


  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 1000,
                           MDES = 0.1, target.power = 0.99,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 100, numCovar.2 = 0,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  calcJ <- pump_sample( design = "d2.2_m2rc",
                        typesample = "J",
                        power.definition = "min1",
                        MTP = "Holm",
                        M = 4,
                        nbar = 1000,
                        MDES = 0.40, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                        rho = 0.2)

  calcJ
  expect_true( !is.na( calcJ$ss.results$`Sample size` ) )


  pp = pump_power( design = "d2.2_m2rc",
                   MTP = "Holm",
                   M = 4,
                   J = calcJ$ss.results$`Sample size` - 1,
                   nbar = 1000,
                   MDES = 0.40,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2,  tnum=1000)
  pp
  expect_true( pp[2,"min1"] <= 0.80 )

  pp = pump_power( design = "d2.2_m2rc",
                   MTP = "Holm",
                   M = 4,
                   J = calcJ$ss.results$`Sample size`,
                   nbar = 1000,
                   MDES = 0.40,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2, tnum=1000 )
  pp
  expect_true( pp[2,"min1"] >= 0.80 )
})


test_that("testing of d3.2_m3rr2rc", {


  traw <- pump_sample_raw( design = "d3.2_m3rr2rc",
                           typesample = "K",
                           nbar = 1000,
                           J = 10,
                           MDES = 0.40, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.05, omega.3 = 0.5 )
  expect_true( !is.na( traw ) )

})

