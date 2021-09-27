
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
    power.definition = "D1indiv",
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



test_that("Bonferroni for non individual power", {

  p <- pump_power(  design = "d2.1_m2fc",
                    MTP = "Bonferroni",
                    J = 10,
                    nbar = 200,
                    M = 3,
                    MDES = rep(0.05, 3),
                    Tbar = 0.50, alpha = 0.05,
                    numCovar.1 = 5, numCovar.2 = 1,
                    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4 )
  p

  ss <- pump_sample(    design = "d2.1_m2fc",
                        MTP = "Bonferroni",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "min1",
                        M = 3,
                        MDES = 0.05, target.power = p$min1[2],
                        tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                        rho = 0.4 )


  expect_equal(ss$`Sample size`, 10, tol = 1)
} )



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
                    J = ss2$`Sample size`,
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

  set.seed( 20303 )
  pow_ref <- pump_power( design = "d2.2_m2rc",
                         MTP = "Holm",
                         M = 4,
                         J = 10,
                         nbar = 10000,
                         MDES = rep( 0.40, 4 ),
                         Tbar = 0.50, alpha = 0.05,
                         numCovar.1 = 5, numCovar.2 = 1,
                         R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                         rho = 0.2, tnum = 10000)

  pow_ref

  # this converges, but not to the correct value
  # because the power curve is too flat
  set.seed( 20303 )
  expect_warning( nbar1 <- pump_sample(
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
  nbar1
  expect_true( !is.na( nbar1$`Sample size` ) )


  # same problem happens with logit
  set.seed( 20303 )
  expect_warning( nbar2 <- pump_sample(
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
    rho = 0.2, just.result.table = FALSE, use.logit = TRUE ) )
  nbar2
  expect_true( !is.na( nbar2$`Sample size` ) )

  # Now an infeasible calculation where the correlation makes min1 not able to
  # achieve power, even though independence would.
  set.seed( 443434344 )
  expect_warning(nbar3 <- pump_sample( design = "d2.2_m2rc",
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
  nbar3
  expect_true( is.na( nbar3$`Sample size` ) )

  # same happens with logit
  set.seed( 443434344 )
  expect_warning(nbar4 <- pump_sample( design = "d2.2_m2rc",
                                       typesample = "nbar",
                                       power.definition = "min1",
                                       MTP = "Holm",
                                       M = 4,
                                       J = 10,
                                       MDES = 0.39, target.power = 0.80, tol = 0.01,
                                       Tbar = 0.50, alpha = 0.05,
                                       numCovar.1 = 5, numCovar.2 = 1,
                                       R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                                       rho = 0.2, max.tnum = 200,
                                       use.logit = TRUE,
                                       just.result.table = FALSE ) )
  nbar4
  expect_true( is.na( nbar4$`Sample size` ) )
})




test_that("Sample with different correlations", {

    # zero correlation
    p <- pump_power(  design = "d2.1_m2fc",
                      MTP = "Holm",
                      J = 10,
                      nbar = 200,
                      M = 20,
                      MDES = rep(0.05, 20),
                      Tbar = 0.50, alpha = 0.05,
                      numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0 )
    p

    ss <- pump_sample(    design = "d2.1_m2fc",
                          MTP = "Holm",
                          typesample = "J",
                          nbar = 200,
                          power.definition = "min1",
                          M = 20,
                          MDES = 0.05, target.power = p$min1[2],
                          tol = 0.01,
                          Tbar = 0.50, alpha = 0.05,
                          numCovar.1 = 5, numCovar.2 = 1,
                          R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                          rho = 0 )


    expect_equal(ss$`Sample size`, 10, tol = 1)


    # high correlation
    p <- pump_power(  design = "d2.1_m2fc",
                      MTP = "Holm",
                      J = 10,
                      nbar = 200,
                      M = 20,
                      MDES = rep(0.05, 20),
                      Tbar = 0.50, alpha = 0.05,
                      numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0.95 )
    p

    ss <- pump_sample(    design = "d2.1_m2fc",
                          MTP = "Holm",
                          typesample = "J",
                          nbar = 200,
                          power.definition = "min1",
                          M = 20,
                          MDES = 0.05, target.power = p$min1[2],
                          tol = 0.01,
                          Tbar = 0.50, alpha = 0.05,
                          numCovar.1 = 5, numCovar.2 = 1,
                          R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
                          rho = 0.95 )


    expect_equal(ss$`Sample size`, 10, tol = 1)

} )

test_that("No adjustment", {

  nbar <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'Bonferroni',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = 0.8,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  )

  nbar <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'None',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = 0.8,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  )

  expect_error(nbar <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'None',
    power.definition = 'complete',
    typesample = 'nbar',
    target.power = 0.8,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  ))

})

