
# Test specific designs for test_pump_sample

# library( pum )
# library( testthat )


test_that("testing of d2.2_m2rc", {

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
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  nbar1

  expect_equal(50, nbar1$`Sample size`, tol = 1)
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

  set.seed( 524235325 )
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
  expect_equal(30, J1$`Sample size`, tol = 5)

  # increasing start.tnum produces results closer to the truth
  set.seed( 524235325 )
  expect_warning( J2 <- pump_sample(
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    start.tnum = 1000) )
  J2
  expect_equal(30, J2$`Sample size`, tol = 1)

  # comparing to logit
  set.seed( 524235325 )
  expect_warning( J3 <- pump_sample(
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    use.logit = TRUE) )
  J3
  expect_equal(30, J3$`Sample size`, tol = 2)

  set.seed( 524235326 )
  expect_warning( J4 <- pump_sample(
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    use.logit = TRUE) )
  J4
  expect_equal(30, J4$`Sample size`, tol = 2)

  # sometimes this converges
  set.seed( 524235325 )
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar1
  expect_true(!is.na(nbar1$ss.results$`Sample size`))
  expect_equal(50, nbar1$ss.results$`Sample size`, tol = 2)

  # sometimes it doesn't (only difference is a new seed)
  set.seed( 524235330 )
  expect_warning(nbar2 <- pump_sample(
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar2
  expect_true(is.na(nbar2$ss.results$`Sample size`))


  # also does not converge with logit
  set.seed( 524235330 )
  expect_warning(nbar3 <- pump_sample(
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    use.logit = TRUE,
    just.result.table = FALSE))
  nbar3
  expect_true(is.na(nbar3$ss.results$`Sample size`))

  # but more iterations fixes it
  set.seed( 524235330 )
  expect_warning(nbar4 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    start.tnum = 1000,
    K = 15,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar4
  expect_true(!is.na(nbar4$ss.results$`Sample size`))
  expect_equal(50, nbar4$ss.results$`Sample size`, tol = 2)

  # more iterations does not fix for logit
  set.seed( 524235330 )
  expect_warning(nbar5 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    start.tnum = 1000,
    K = 15,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE,
    use.logit = TRUE))
  nbar5
  expect_true(!is.na(nbar5$ss.results$`Sample size`))
  expect_equal(50, nbar5$ss.results$`Sample size`, tol = 2)

  # decreasing the maximum sample size fixes it
  set.seed( 524235330 )
  expect_warning(nbar6 <- pump_sample(
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
  nbar6
  expect_true(!is.na(nbar6$ss.results$`Sample size`))
  expect_equal(50, nbar6$ss.results$`Sample size`, tol = 10)

  # decreasing max sample size for logit makes it converge
  # to the wrong value!
  set.seed( 524235330 )
  expect_warning(nbar7 <- pump_sample(
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE,
    use.logit = TRUE))
  nbar7
  expect_true(!is.na(nbar7$ss.results$`Sample size`))
  expect_equal(50, nbar7$ss.results$`Sample size`, tol = 30)

  # decreasing max sample size and increasing start.tnum does not help
  set.seed( 524235330 )
  expect_warning(nbar8 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    max_sample_size_nbar = 1000,
    start.tnum = 1000,
    K = 15,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE,
    use.logit = TRUE))
  nbar8
  expect_true(!is.na(nbar8$ss.results$`Sample size`))
  expect_equal(50, nbar8$ss.results$`Sample size`, tol = 30)
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

  # sometimes this doesn't converge
  set.seed( 245444 )
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar1
  expect_true(is.na(nbar1$ss.results$`Sample size`))

  # logit doesn't converge either
  set.seed( 245444 )
  expect_warning(nbar2 <- pump_sample(
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE,
    use.logit = TRUE))
  nbar2
  expect_true(is.na(nbar2$ss.results$`Sample size`))

  # decreasing max sample size helps, but it is
  # still pretty far from the true value!
  set.seed( 245444 )
  expect_warning(nbar3 <- pump_sample(
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar3
  expect_equal(50, nbar3$ss.results$`Sample size`, tol = 15)

  # logit has the same pattern
  set.seed( 245444 )
  expect_warning(nbar4 <- pump_sample(
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE,
    use.logit = TRUE))
  nbar4
  expect_equal(50, nbar4$ss.results$`Sample size`, tol = 20)

  # increasing start.tnum is not enough
  set.seed( 245444 )
  expect_warning(nbar5 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    start.tnum = 1000,
    K = 10,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar5
  expect_true(is.na(nbar5$ss.results$`Sample size`))


  # trying both doesn't improve the closeness to the truth
  set.seed( 245444 )
  expect_warning(nbar6 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    max_sample_size_nbar = 1000,
    start.tnum = 1000,
    K = 10,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5, just.result.table = FALSE))
  nbar6
  expect_equal(50, nbar6$ss.results$`Sample size`, tol = 12)


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

  # does not converge
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

  expect_true(is.na(J1$`Sample size`))

  # does not converge with logit
  J2 <- pump_sample(
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    use.logit = TRUE)
  expect_true(is.na(J2$`Sample size`))

  # does not converge
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

  expect_true(is.na(nbar1$`Sample size`))
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
