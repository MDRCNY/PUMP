
# Test specific designs for test_pump_sample

# library( pum )
# library( testthat )

#------------------------------------------------------------------#
# test pump sample raw
#------------------------------------------------------------------#

test_that("testing of d2.2_m2rc raw", {

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
  expect_true( !is.na( calcJ$`Sample size` ) )


  pp = pump_power( design = "d2.2_m2rc",
                   MTP = "Holm",
                   M = 4,
                   J = calcJ$`Sample size` - 1,
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
                   J = calcJ$`Sample size`,
                   nbar = 1000,
                   MDES = rep( 0.40, 4),
                   Tbar = 0.50, alpha = 0.05,
                   numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2, tnum=1000 )
  pp
  expect_true( pp[2,"min1"] >= 0.80 )
})

test_that("testing of d3.2_m3rr2rc raw", {

  traw <- pump_sample_raw( design = "d3.2_m3rr2rc",
                           typesample = "K",
                           nbar = 1000,
                           J = 10,
                           MDES = 0.40, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = 0.05, ICC.3 = 0.05,
                           omega.3 = 0.5 )
  traw
  expect_true( !is.na( traw$K ) )
  expect_true( !is.na( traw$df ) )

})

#------------------------------------------------------------------
# two level models
#------------------------------------------------------------------

#----------------------------------------
# d2.1_m2fc
#----------------------------------------

test_that("testing of d2.1_m2fc", {

  set.seed(8598)

  pp1 <- pump_power(
    design = "d2.1_m2fc",
    MTP = 'Holm',
    nbar = 50,
    J = 60,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)

  J1 <- pump_sample(
    design = "d2.1_m2fc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp1$D1indiv[2],
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample size`, 60, tolerance = 0.1)

  # converges
  set.seed(8598)
  nbar1 <- pump_sample(
    design = "d2.1_m2fc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp1$D1indiv[2],
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
    rho = 0.2)
  nbar1
  expect_equal(50, nbar1$`Sample size`, tol = 0.1)

})

#----------------------------------------
# d2.1_m2ff
#----------------------------------------

test_that("testing of d2.1_m2ff", {

  set.seed(8598)

  pp1 <- pump_power(
    design = "d2.1_m2ff",
    MTP = 'Holm',
    nbar = 50,
    J = 60,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)

  J1 <- pump_sample(
    design = "d2.1_m2ff",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp1$D1indiv[2],
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample size`, 60, tol = 0.1)

  set.seed(8598)
  nbar1 <- pump_sample(
    design = "d2.1_m2ff",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp1$D1indiv[2],
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  nbar1
  expect_equal(nbar1$`Sample size`, 50, tol = 0.1)

})

#----------------------------------------
# d2.1_m2fr
#----------------------------------------

test_that("testing of d2.1_m2fr", {

  set.seed(8598)

  pp1 <- pump_power(
    design = "d2.1_m2fr",
    MTP = 'Holm',
    nbar = 50,
    J = 60,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)

  J1 <- pump_sample(
    design = "d2.1_m2ff",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp1$D1indiv[2],
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample size`, 60, tol = 0.1)

  nbar1 <- pump_sample(
    design = "d2.1_m2fr",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp1$D1indiv[2],
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  nbar1
  expect_equal(nbar1$`Sample size`, 50, tol = 0.1)

})

#----------------------------------------
# d2.2_m2rc
#----------------------------------------

test_that("testing of d2.2_m2rc", {

  set.seed(8598)

  pp1 <- pump_power(
    design = "d2.2_m2rc",
    MTP = 'Holm',
    nbar = 50,
    J = 60,
    M = 3,
    MDES = rep(0.125, 3),
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)

  J1 <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp1$D1indiv[2],
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample size`, 60, tol = 0.1)

  nbar1 <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp1$D1indiv[2],
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  nbar1
  expect_equal(nbar1$`Sample size`, 50, tol = 0.1)

})

#------------------------------------------------------------------
# three level models
#------------------------------------------------------------------

#----------------------------------------
# d3.1_m3rr2rr
#----------------------------------------

test_that("testing of d3.1_m3rr2rr", {

  set.seed( 524235326 )
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

  set.seed( 524235326 )
  K1 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'K',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    J = 30,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5)
  K1
  expect_equal(K1$`Sample size`, 15, tol = 0.1)

  # converges
  set.seed( 524235326 )
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
    omega.2 = 0.1, omega.3 = 0.1,
    rho = 0.5,
    just.result.table = FALSE ) )
  J1
  expect_equal(J1$`Sample size`, 30, tol = 0.2)

  # in this case it is very flat
  set.seed( 524235327 )
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
    just.result.table = FALSE) )
  J2
  expect_equal(J2$`Sample size`, 30, tol = 0.3)

  # decreasing tolerance and increasing max steps does help!
  set.seed( 524235327 )
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
    just.result.table = FALSE,
    tol = 0.005, max.steps = 40, start.tnum = 2000, max.tnum = 4000) )
  J3
  expect_equal(J3$`Sample size`, 30, tol = 0.2)

  # decreasing max sample size does not help
  set.seed( 524235327 )
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
    just.result.table = FALSE,
    max_sample_size_JK = 100) )
  J4
  expect_equal(J3$`Sample size`, 30, tol = 0.4)

  # sometimes this converges
  set.seed( 524235327 )
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar1
  expect_true(!is.na(nbar1$`Sample size`))
  # but not very close!
  expect_equal(50, nbar1$`Sample size`, tol = 0.4)

  # sometimes it doesn't (only difference is a new seed)
  set.seed( 524235325 )
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
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar2
  expect_true(is.na(nbar2$`Sample size`))

  # decreasing the maximum sample size alone does not fix it
  set.seed( 524235325 )
  expect_warning(nbar5 <- pump_sample(
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
    max_sample_size_nbar = 500,
    just.result.table = FALSE))
  nbar5
  expect_true(is.na(nbar5$`Sample size`))

  # it converges with more iterations
  # but it is very flat!
  set.seed( 524235325 )
  expect_warning(nbar4 <- pump_sample(
    design = "d3.1_m3rr2rr",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    start.tnum = 10000,
    K = 15,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0.1, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar4
  expect_true(!is.na(nbar4$`Sample size`))
  expect_equal(nbar4$`Sample size`, 50, tol = 0.5)

  # decreasing max sample size helps get a smaller estiamte
  set.seed( 524235325 )
  expect_warning(nbar5 <- pump_sample(
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
    max_sample_size_nbar = 100,
    just.result.table = FALSE))
  nbar5
  expect_true(!is.na(nbar5$`Sample size`))
  expect_equal(nbar5$`Sample size`, 50, tol = 0.3)

  # decreasing max sample size helps get a smaller estiamte
  # and decreasing tolerance
  set.seed( 524235325 )
  expect_warning(nbar6 <- pump_sample(
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
    max_sample_size_nbar = 100,
    tol = 0.005,
    just.result.table = FALSE))
  nbar6
  expect_true(!is.na(nbar6$`Sample size`))
  expect_equal(nbar6$`Sample size`, 50, tol = 0.2)


  # if we go below the true value, does not converge!
  set.seed( 524235325 )
  expect_warning(nbar7 <- pump_sample(
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
    max_sample_size_nbar = 40,
    just.result.table = FALSE))
  nbar7
  expect_true(is.na(nbar7$`Sample size`))
})


#----------------------------------------
# d3.2_m3ff2rc
#----------------------------------------

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

  K1 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'K',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    nbar = 50,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)
  K1
  expect_equal(K1$`Sample size`, 10, tol = 0.1)

  J1 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    nbar = 50,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)
  J1
  expect_equal(J1$`Sample size`, 30, tol = 0.1)

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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar1
  expect_true(is.na(nbar1$`Sample size`))

  # decreasing max sample size helps, but it is
  # very flat
  set.seed( 245444 )
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar2
  expect_equal(nbar2$`Sample size`, 50, tol = 0.1)

})


#----------------------------------------
# d3.2_m3ff2rc
#----------------------------------------

test_that("testing of d3.2_m3ff2rc", {
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

  set.seed( 245444 )
  K1 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'K',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    nbar = 50,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)
  K1
  expect_equal(K1$`Sample size`, 10, tol = 0.1)

  # converges
  set.seed( 245444 )
  J1 <- pump_sample(
    design = "d3.2_m3ff2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    nbar = 50,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)
  J1
  expect_equal(J1$`Sample size`, 30, tol = 0.1)

  # does not converge
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
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar1
  expect_true(is.na(nbar1$`Sample size`))

  # does converge with smaller max sample size
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
    max_sample_size_nbar = 1000,
    just.result.table = FALSE))
  nbar2
  expect_equal(nbar2$`Sample size`, 50, tol = 0.1)

})


#----------------------------------------
# d3.2_m3rr2rc
#----------------------------------------

test_that("testing of d3.2_m3rr2rc", {
  set.seed( 245444 )

  pp1 <- pump_power(
    design = "d3.2_m3rr2rc",
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

  set.seed( 245444 )
  K1 <- pump_sample(
    design = "d3.2_m3rr2rc",
    typesample = 'K',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    nbar = 50,
    J = 30,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)
  K1
  expect_equal(K1$`Sample size`, 10, tol = 0.1)

  # this converges
  set.seed( 245444 )
  J1 <- pump_sample(
    design = "d3.2_m3rr2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    nbar = 50,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5)
  J1
  expect_equal(J1$`Sample size`, 30, tol = 0.1)

  # this converges to a value that is way too high!
  set.seed( 245444 )
  nbar1 <- expect_warning(pump_sample(
    design = "d3.2_m3rr2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    J = 30,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    just.result.table = FALSE))
  nbar1
  expect_true(!is.na(nbar1$`Sample size`))

  # this converges to a value that is less high
  # but still very high
  set.seed( 245444 )
  nbar2 <- expect_warning(pump_sample(
    design = "d3.2_m3rr2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    J = 30,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    max_sample_size_nbar = 1000,
    just.result.table = FALSE))
  nbar2
  expect_true(!is.na(nbar2$`Sample size`))

  # decreasing tolerance makes it not converge
  set.seed( 245444 )
  nbar3 <- expect_warning(pump_sample(
    design = "d3.2_m3rr2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    J = 30,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    tol = 0.005, max.steps = 40,
    just.result.table = FALSE))
  nbar3
  expect_true(is.na(nbar3$`Sample size`))


  # decreasing tolerance and decreasing max sample size
  # still gives a big over estimate
  set.seed( 245444 )
  nbar4 <- expect_warning(pump_sample(
    design = "d3.2_m3rr2rc",
    typesample = 'nbar',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    J = 30,
    K = 10,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.1,
    ICC.2 = 0.2, ICC.3 = 0.2,
    omega.2 = 0, omega.3 = 0.1, rho = 0.5,
    tol = 0.005, max.steps = 60, max.tnum = 4000, start.tnum = 2000,
    max_sample_size_nbar = 1000,
    just.result.table = FALSE))
  nbar4
  expect_true(!is.na(nbar4$`Sample size`))

  # try decreasing max sample size more
  set.seed( 245444 )
  nbar5 <- expect_warning(pump_sample(
    design = "d3.2_m3rr2rc",
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
    tol = 0.005, max.steps = 60, max.tnum = 4000, start.tnum = 2000,
    max_sample_size_nbar = 200,
    just.result.table = FALSE))
  nbar5
  expect_true(!is.na(nbar5$`Sample size`))

  # decreasing max sample size more
  # produces the best answer
  set.seed( 245444 )
  nbar6 <- expect_warning(pump_sample(
    design = "d3.2_m3rr2rc",
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
    tol = 0.005, max.steps = 60, max.tnum = 4000, start.tnum = 2000,
    max_sample_size_nbar = 75,
    just.result.table = FALSE))
  nbar6
  expect_equal(nbar6$`Sample size`, 50, tol = 0.1)

})

#----------------------------------------
# d3.3_m3rc2rc
#----------------------------------------
# does not converge for nbar

test_that("testing of d3.3_m3rc2rc", {

  set.seed(2344)

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


  K1 <- pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'K',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    J = 40,
    nbar = 50,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5)
  K1
  expect_equal(K1$`Sample size`, 20, tol = 0.1)

  # converges but is very flat
  set.seed( 245444 )
  J1 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 20,
    nbar = 50,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    just.result.table = FALSE))
  J1
  expect_true(!is.na(J1$`Sample size`))


  # getting closer!
  set.seed( 245444 )
  J2 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 20,
    nbar = 50,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    tol = 0.005, max_sample_size_JK = 80,
    just.result.table = FALSE))
  J2
  expect_true(!is.na(J2$`Sample size`))

  # let's get closer
  set.seed( 245444 )
  J3 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 20,
    nbar = 50,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    start.tnum = 2000, max.tnum = 4000,
    tol = 0.005, max_sample_size_JK = 80,
    just.result.table = FALSE))
  J3
  expect_true(!is.na(J3$`Sample size`))

  # let's get there!
  set.seed( 245444 )
  J4 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    typesample = 'J',
    MTP = 'Holm',
    target.power = pp1$D1indiv[2],
    power.definition = 'D1indiv',
    K = 20,
    nbar = 50,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    start.tnum = 2000, max.tnum = 4000,
    tol = 0.005, max_sample_size_JK = 50,
    just.result.table = FALSE))
  J4
  expect_true(!is.na(J4$`Sample size`))
  expect_equal(J4$`Sample size`, 40, tol = 0.2)

  # does not converge
  set.seed( 245444 )
  nbar1 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    power.definition = 'D1indiv',
    target.power = pp1$D1indiv[2],
    typesample = 'nbar',
    MTP = 'Holm',
    K = 20,
    J = 40,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    just.result.table = FALSE))
  nbar1
  expect_true(is.na(nbar1$`Sample size`))

  # does not converge
  set.seed( 245444 )
  nbar2 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    power.definition = 'D1indiv',
    target.power = pp1$D1indiv[2],
    typesample = 'nbar',
    MTP = 'Holm',
    K = 20,
    J = 40,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    max_sample_size_nbar = 1000,
    start.tnum = 5000, max.tnum = 8000,
    just.result.table = FALSE))
  nbar2
  expect_true(is.na(nbar2$`Sample size`))

  # does not converge
  set.seed( 245444 )
  nbar3 <- expect_warning(pump_sample(
    design = "d3.3_m3rc2rc",
    power.definition = 'D1indiv',
    target.power = pp1$D1indiv[2],
    typesample = 'nbar',
    MTP = 'Holm',
    K = 20,
    J = 40,
    M = 3,
    MDES = 0.25,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1, numCovar.3 = 1,
    R2.1 = 0.1, R2.2 = 0.1, R2.3 = 0.1,
    ICC.2 = 0.1, ICC.3 = 0.1,
    omega.2 = 0, omega.3 = 0, rho = 0.5,
    max_sample_size_nbar = 200,
    start.tnum = 5000, max.tnum = 8000,
    just.result.table = FALSE))
  nbar3
  expect_true(is.na(nbar3$`Sample size`))
})


#----------------------------------------
# lower limit
#----------------------------------------

test_that( "testing of lower limit", {
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
  expect_true( pp$`min1 power` > 0.50 )
} )
