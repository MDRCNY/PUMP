
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
    Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
    numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  )

  calcnbar <- pump_sample_raw( design = "d2.2_m2rc",
                               typesample = "nbar",
                               J = 10,
                               MDES = 0.05, target.power = 0.80,
                               Tbar = 0.50, alpha = 0.05,
                               numCovar.1 = 5, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcnbar
  expect_true( is.na( calcnbar$nbar ) )


  calcJ <- pump_sample_raw( design = "d2.1_m2fc",
                            typesample = "J",
                            nbar = 258,
                            MDES = 0.05, target.power = 0.80,
                            Tbar = 0.50, alpha = 0.05,
                            numCovar.1 = 5, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcJ
  expect_true( !is.na( calcJ$J ) )

  calcn <- pump_sample_raw( design = "d2.1_m2fc",
                            typesample = "nbar",
                            J = calcJ$J,
                            MDES = 0.05, target.power = 0.80,
                            Tbar = 0.50, alpha = 0.05,
                            numCovar.1 = 5, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcn
  expect_true( abs( calcn$nbar - 258 ) < 2 )

  calcJ2 <- pump_sample_raw( design = "d2.1_m2fc",
                             typesample = "J",
                             nbar = calcn$nbar,
                             MDES = 0.05, target.power = 0.80,
                             Tbar = 0.50, alpha = 0.05,
                             numCovar.1 = 5, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcJ2
  expect_true( abs( calcJ$J - calcJ2$J ) < 1 )

  calcn2 <- pump_sample_raw( design = "d2.1_m2fc",
                             typesample = "nbar",
                             J = calcJ2$J,
                             MDES = 0.05, target.power = 0.80,
                             Tbar = 0.50, alpha = 0.05,
                             numCovar.1 = 5, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcn2
  expect_true( abs( calcn$nbar - calcn2$nbar ) < 2 )
})


test_that( "optimize_power solves", {

  set.seed( 3042424 )
  op_pow <- pum:::optimize_power(
    MTP = "Holm", nbar=200,
    power.definition="D1indiv",
    design = "d2.1_m2fc", search.type = "J",
    start.low = 56, start.high = 75,
    start.tnum = 200,
    M = 3,
    MDES = 0.05, target.power = 0.80, tol = 0.01,
    Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
    rho = 0.4, max.tnum = 400, final.tnum = 2000
  )
  op_pow
  expect_true( ncol( op_pow ) == 6 )
  expect_true( all( op_pow$w <= 2000 ) )
  expect_true( max( op_pow$w ) == 2000 )

})


test_that("pump_sample 2 level/2 level", {
  ss2 <- pump_sample(   design = "d2.1_m2fc",
                        MTP = "Holm",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                        rho = 0.4, just.result.table = FALSE )
  ss2

  p2 <- pump_power( design = "d2.1_m2fc",
                    MTP = "Holm",
                    J = ss2$ss.results$`Sample size`,
                    nbar = 200,
                    M = 3,
                    MDES = rep(0.05, 3),
                    Tbar = 0.50, alpha = 0.05,
                    numCovar.1 = 5, numCovar.2 = 1,
                    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4 )

  p2
  expect_equal( p2[ 2, "indiv.mean" ], 0.80, tol = 0.02 )
} )


test_that("sample search when one end is missing", {

  pow_ref <- pump_power( design = "d2.2_m2rc",
                         MTP = "Holm",
                         M = 4,
                         J = 10,
                         nbar = 10000,
                         MDES = rep( 0.40, 4 ),
                         Tbar = 0.50, alpha = 0.05,
                         numCovar.1 = 5, numCovar.2 = 1,
                         R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                         rho = 0.2, tnum = 200)

  pow_ref

  set.seed( 20303 )
  expect_warning( calcnbar <- pump_sample(
    design = "d2.2_m2rc",
    typesample = "nbar",
    power.definition = "min1",
    MTP = "Holm",
    M = 4,
    J = 10,
    MDES = 0.40, target.power = pow_ref$min1[2], tol = 0.01,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
    rho = 0.2, just.result.table = FALSE ) )

  calcnbar
  expect_true( !is.na( calcnbar$ss.results$`Sample size` ) )

  # Now an infeasible calculation where the correlation makes min1 not able to
  # achieve power, even though independence would.
  set.seed( 443434344 )
  expect_warning(calcnbar <- pump_sample( design = "d2.2_m2rc",
                                          typesample = "nbar",
                                          power.definition = "min1",
                                          MTP = "Holm",
                                          M = 4,
                                          J = 10,
                                          MDES = 0.39, target.power = 0.80, tol = 0.01,
                                          Tbar = 0.50, alpha = 0.05,
                                          numCovar.1 = 5, numCovar.2 = 1,
                                          R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                                          rho = 0.2, max.tnum = 200, just.result.table = FALSE ) )
  calcnbar
  expect_true( is.na( calcnbar$ss.results$`Sample size` ) )
})






test_that("further testing of d2.2_m2rc", {

  set.seed( 101010 )

  traw <- pump_sample_raw(
    design = "d2.2_m2rc",
    typesample = "J",
    nbar = 1000,
    MDES = 0.40, target.power = 0.80,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 10,
                           MDES = 0.40, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 10,
                           MDES = 0.01, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 1000,
                           MDES = 0.001, target.power = 0.99,
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw


  traw <- pump_sample_raw( design = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 1000,
                           MDES = 0.1, target.power = 0.99,
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 100, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  calcJ <- pump_sample( design = "d2.2_m2rc",
                        typesample = "J",
                        power.definition = "min1",
                        MTP = "Holm",
                        M = 4,
                        nbar = 1000,
                        MDES = 0.40, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                        max.tnum = 1000,
                        rho = 0.2, just.result.table = FALSE)

  calcJ
  expect_true( !is.na( calcJ$ss.results$`Sample size` ) )


  pp = pump_power( design = "d2.2_m2rc",
                   MTP = "Holm",
                   M = 4,
                   J = calcJ$ss.results$`Sample size` - 1,
                   nbar = 1000,
                   MDES = rep(0.40, 4),
                   Tbar = 0.50, alpha = 0.05,
                   numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2,  tnum=1000)
  pp
  expect_true( pp[2,"min1"] <= 0.80 )

  pp = pump_power( design = "d2.2_m2rc",
                   MTP = "Holm",
                   M = 4,
                   J = calcJ$ss.results$`Sample size`,
                   nbar = 1000,
                   MDES = rep( 0.40, 4),
                   Tbar = 0.50, alpha = 0.05,
                   numCovar.1 = 5, numCovar.2 = 1,
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
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.05, omega.3 = 0.5 )
  traw
  expect_true( !is.na( traw$K ) )
  expect_true( !is.na( traw$df ) )

})


test_that("testing of d2.2_m2rc", {

  pp1 <- pump_power(
    design = "d2.2_m2rc",
    MTP = 'Bonferroni',
    nbar = 50,
    J = 60,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)

  nbar1 <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'Bonferroni',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp1$D1indiv[2],
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2, just.result.table = FALSE)
  nbar1

  expect_equal(50, nbar1$`Sample size`[[2]], tol = 1)
})

test_that("testing of d3.1_m3rr2rr", {

  pp1 <- pump_power(
    design = "d3.1_m3rr2rr",
    MTP = 'Holm',
    nbar = 50,
    K = 15,
    J = 30,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5)

  expect_warning( J1 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 15,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5) )
  J1
  expect_equal(30, J1$`Sample size`, tol = 3)

  # expect_true(is.na(J1$ss.results$`Sample size`))
  #
  # expect_warning(J2 <- pump_sample(
  #     design = "d3.1_m3rr2rr",
  #     typesample = 'J',
  #     MTP = 'Holm',
  #     target.power = pp1$D1indiv[2],
  #     power.definition = 'D1indiv',
  #     max_sample_size_JK = 100,
  #     K = 15,
  #     nbar = 50,
  #     M = 3,
  #     MDES = 0.125,
  #     Tbar = 0.5, alpha = 0.05,
  #     numCovar.1 = 1, numCovar.2 = 1,
  #     R2.1 = 0.1, R2.2 = 0.1,
  #     ICC.2 = 0.2, ICC.3 = 0.2,
  #     omega.2 = 0.1, omega.3 = 0.1, rho = 0.5))
  #
  # expect_true(!is.na(J2$ss.results$`Sample size`))
  # expect_equal(30, J2$ss.results$`Sample size`, tol = 5)

  expect_warning(nbar1 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 15,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5))
  nbar1
  expect_true(is.na(nbar1$`Sample size`))

  set.seed( 524235325 )
  expect_warning(nbar2 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    max_sample_size_nbar = 1000,
    K = 15,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar2
  expect_true(!is.na(nbar2$ss.results$`Sample size`))
  expect_equal(50, nbar2$ss.results$`Sample size`, tol = 2)
})


test_that("testing of d3.1_m3ff2rr", {
  set.seed( 245444 )

  pp1 <- pump_power(
    design = "d3.2_m3ff2rc",
    MTP = 'Holm',
    nbar = 50,
    K = 10,
    J = 30,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)

  expect_warning(nbar1 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 10,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5))
  nbar1
  expect_true(is.na(nbar1$`Sample size`))

  expect_warning(nbar2 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    max_sample_size_nbar = 1000,
    K = 10,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5))

  # expect_equal(50, nbar1$ss.results$`Sample size`, tol = 5)
})


test_that("testing of d3.3_m3rc2rc", {

  pp1 <- pump_power(
    design = "d3.3_m3rc2rc",
    MTP = 'Holm',
    nbar = 50,
    K = 20,
    J = 40,
    M = 3,
    MDES = rep(0.25, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5)

  J1 <- pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 10,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)

  expect_true(is.na(J1$ss.results$`Sample size`))

  nbar1 <- pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 10,
    J = 40,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)

  expect_true(is.na(nbar1$ss.results$`Sample size`))
})




test_that( "testing pump_sample for d3.2_m3ff2rc", {
  # This should hit lower limit (too powerful, want J < 3).
  set.seed( 24553453 )
  expect_warning( pp <- pump_sample(    design = "d3.2_m3ff2rc",
                                        typesample = "J",
                                        MTP = "Holm",
                                        MDES = 0.12,
                                        target.power = 0.50,
                                        power.definition = "min1",
                                        tol = 0.01,
                                        M = 5,
                                        K = 7, # number RA blocks
                                        nbar = 58,
                                        Tbar = 0.50, # prop Tx
                                        alpha = 0.15, # significance level
                                        numCovar.1 = 1, numCovar.2 = 1,
                                        R2.1 = 0.1, R2.2 = 0.7,
                                        ICC.2 = 0.05, ICC.3 = 0.9,
                                        rho = 0.4, # how correlated outcomes are
                                        max.tnum = 200, just.result.table = FALSE ) )
  pp
  expect_true( !is.null( pp ) )
  expect_true( pp$ss.results$`min1 power` > 0.50 )
} )
