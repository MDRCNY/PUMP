
# Test specific designs for test_pump_sample

# library( pum )
# library( testthat )

source( "testing_code.R" )



# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -----------------------------------------------#
# test pump sample raw
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -----------------------------------------------#

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
                        rho = 0.2)

  calcJ
  expect_true( !is.na( calcJ$`Sample.size` ) )


  pp = pump_power( design = "d2.2_m2rc",
                   MTP = "Holm",
                   M = 4,
                   J = calcJ$`Sample.size` - 1,
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
                   J = calcJ$`Sample.size`,
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
  expect_true( !is.na( traw$ss ) )
  expect_true( !is.na( traw$df ) )

})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# --------- two level models --------

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# -------- d2.1_m2fc --------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fc", {

  if ( FALSE ) {
    set.seed(8598)

    pp1 <- pump_power(
      design = "d2.1_m2fc",
      MTP = 'Holm',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05,
      numCovar.1 = 1, numCovar.2 = 1, tnum = 100000,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
    pp1
    pp_power = pp1$D1indiv[[2]]

  }
  pp_power = 0.95162

  J1 <- pump_sample(
    design = "d2.1_m2fc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp_power,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample.size`, 60, tolerance = 0.1)

  # converges
  set.seed(8598)
  nbar1 <- pump_sample(
    design = "d2.1_m2fc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power =pp_power,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
    rho = 0.2)
  nbar1
  expect_equal(50, nbar1$`Sample.size`, tol = 0.1)

})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# -------------  d2.1_m2ff -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2ff", {

  if ( FALSE ) {

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
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2, tnum = 100000)
  pp1
  }

  pp_power = 0.95152
  vals = test_sample_triad(pp_power, 50, 60, NULL, 24322323,
    design = "d2.1_m2ff",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)

  expect_equal(60, vals$J, tol=0.1 )
  expect_equal(50, vals$nbar, tol=0.1 )
  expect_equal( warning_pattern(vals), c(FALSE,FALSE) )


})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d2.1_m2fr -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fr", {

  if ( FALSE ) {
    set.seed(8598)

    pp1 <- pump_power(
      design = "d2.1_m2fr",
      MTP = 'Holm',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05,
      numCovar.1 = 1, numCovar.2 = 1, tnum = 100000,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
    pp1$D1indiv[[2]]
  }
  pp_power = 0.9424

  vals = test_sample_triad(pp_power, nbar=50, J=60, K=NULL, seed=22422422,
                           design = "d2.1_m2ff",
                           MTP = 'Holm',
                           power.definition = 'D1indiv',
                           M = 3,
                           MDES = 0.125,
                           Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  vals[1:2]
  warning_pattern(vals)

  expect_equal(60, vals$J, tol=0.1)
  expect_equal(50, vals$nbar, tol=0.1)
  expect_equal( warning_pattern(vals), c(FALSE,FALSE) )

})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -

# ------   d2.2_m2rc   -------

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -

test_that("testing of d2.2_m2rc", {

  if ( FALSE ) {

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

    pp1
    pp1$D1indiv[2]
  }
  pp_power = 0.688

  J1 <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power =pp_power,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample.size`, 60, tol = 0.1)

  nbar1 <- pump_sample(
    design = "d2.2_m2rc",
    MTP = 'Holm',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power =pp_power,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
  nbar1
  expect_equal(nbar1$`Sample.size`, 50, tol = 0.1)

})
