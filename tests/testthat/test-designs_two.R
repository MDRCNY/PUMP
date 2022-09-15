# library( PUMP )
# library( testthat )

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# --------- two level models --------

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -----------------------------------------------#
# test pump sample raw
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -----------------------------------------------#

skip_on_cran()

test_that("testing of d2.2_m2rc raw", {

  set.seed( 101010 )

  traw <- pump_sample_raw(
    d_m = "d2.2_m2rc",
    typesample = "J",
    nbar = 1000,
    MDES = 0.40, target.power = 0.80,
    Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 5, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  traw <- pump_sample_raw( d_m = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 10,
                           MDES = 0.40, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  traw <- pump_sample_raw( d_m = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 10,
                           MDES = 0.01, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  traw <- pump_sample_raw( d_m = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 1000,
                           MDES = 0.001, target.power = 0.99,
                           Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw


  traw <- pump_sample_raw( d_m = "d2.2_m2rc",
                           typesample = "J",
                           nbar = 1000,
                           MDES = 0.1, target.power = 0.99,
                           Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                           numCovar.1 = 100, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  traw

  set.seed( 1041010 )
  
  calcJ <- pump_sample( d_m = "d2.2_m2rc",
                        typesample = "J",
                        power.definition = "min1",
                        MTP = "HO",
                        M = 4,
                        nbar = 1000,
                        MDES = 0.40, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                        numCovar.1 = 5, numCovar.2 = 1,
                        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                        tnum = 1000,
                        rho = 0.2)

  calcJ
  expect_true( !is.na( calcJ$`Sample.size` ) )

  pp <- pump_power( d_m = "d2.2_m2rc",
                   MTP = "HO",
                   M = 4,
                   J = calcJ$`Sample.size` - 1,
                   nbar = 1000,
                   MDES = rep(0.40, 4),
                   Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                   numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2,  tnum=1000)
  pp
  expect_true( pp[2,"min1"] <= 0.80 )

  pp <- pump_power( d_m = "d2.2_m2rc",
                   MTP = "HO",
                   M = 4,
                   J = calcJ$`Sample.size`,
                   nbar = 1000,
                   MDES = rep( 0.40, 4),
                   Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                   numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2, tnum=1000 )
  pp
  expect_true( pp[2,"min1"] >= 0.80 )
})

test_that("testing of d3.2_m3rr2rc raw", {

  traw <- pump_sample_raw( d_m = "d3.2_m3rr2rc",
                           typesample = "K",
                           nbar = 1000,
                           J = 10,
                           MDES = 0.40, target.power = 0.80,
                           Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7,
                           ICC.2 = 0.05, ICC.3 = 0.05,
                           omega.3 = 0.5 )
  traw
  expect_true( !is.na( traw$ss ) )
  expect_true( !is.na( traw$df ) )

})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# -------- d1.1_m1c --------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d1.1_m1c", {
    
    pp <- pump_sample( d_m = "d1.1_m1c",
                       MTP = c("BF"),
                       MDES = 0.125,
                       power.definition = 'D1indiv',
                       target.power = 0.8,
                       typesample = 'nbar',
                       M = 5,
                       Tbar = 0.50, # prop Tx
                       alpha = 0.05, # significance level
                       numCovar.1 = 5,
                       R2.1 = 0.1, 
                       rho = 0.4, # how correlated outcomes are
                       tnum = 1000
    )
    
  expect_true(!is.null(pp))
})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# -------- d2.1_m2fc --------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fc one-tailed", {

  if ( FALSE ) {
    set.seed(8598)

    pp1 <- pump_power(
      d_m = "d2.1_m2fc",
      MTP = 'HO',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
      numCovar.1 = 1, numCovar.2 = 1, tnum = 100000,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
    pp1
    pp_power <- pp1$D1indiv[2]

  }

  pp_power <- 0.97633

  J1 <- pump_sample(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp_power,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(J1$`Sample.size`, 60, tolerance = 0.1)

  # converges
  set.seed(8598)
  nbar1 <- pump_sample(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp_power,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1, 
    R2.1 = 0.1, ICC.2 = 0.05,
    rho = 0.2)
  nbar1
  expect_equal(50, nbar1$`Sample.size`, tol = 0.1)

  set.seed( 44304044 )
  mdes1 <- pump_mdes(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    power.definition = 'D1indiv',
    target.power = pp_power,
    J = 60,
    nbar = 50,
    M = 3,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05,
    rho = 0.2)
  expect_equal(0.125, mdes1$Adjusted.MDES, tolerance = 0.1)

})


test_that("testing of d2.1_m2fc two-tailed", {

  if ( FALSE ) {
    set.seed(8598)

    pp1 <- pump_power(
      d_m = "d2.1_m2fc",
      MTP = 'HO',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1, numCovar.2 = 1, tnum = 100000,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2)
    pp1
    pp_power <- pp1$D1indiv[2]

  }

  pp_power <- 0.95162

  J1 <- pump_sample(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    power.definition = 'D1indiv',
    typesample = 'J',
    target.power = pp_power,
    nbar = 50,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2)
  J1
  expect_equal(60, J1$`Sample.size`, tolerance = 0.1)

  # converges
  set.seed(8598)
  nbar1 <- pump_sample(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = pp_power,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05,
    rho = 0.2)
  nbar1
  expect_equal(50, nbar1$`Sample.size`, tol = 0.1)

  mdes1 <- pump_mdes(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    power.definition = 'D1indiv',
    target.power = pp_power,
    J = 60,
    nbar = 50,
    M = 3,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05,
    rho = 0.2)
  expect_equal(0.125, mdes1$Adjusted.MDES, tolerance = 0.1)

})




# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# -------------  d2.1_m2ff -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2ff one-tailed", {

  if ( FALSE ) {

    set.seed(8598)

    pp1 <- pump_power(
      d_m = "d2.1_m2ff",
      MTP = 'HO',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2,
      tnum = 100000 )
    pp1
    pp_power <- pp1$D1indiv[2]

  }

  pp_power <- 0.97644

  vals <- test_sample_triad(pp_power, nbar = 50, J = 60, NULL, 24322323,
    d_m = "d2.1_m2ff",
    MTP = 'HO',
    power.definition = 'D1indiv',
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2)

  expect_equal(60, vals$J, tol = 0.1 )
  expect_equal(50, vals$nbar, tol = 0.1 )
  expect_equal( warning_pattern(vals), c(FALSE, FALSE) )
})

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d2.1_m2fr -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fr one-tailed", {

  if ( FALSE ) {
    set.seed(8598)

    pp1 <- pump_power(
      d_m = "d2.1_m2fr",
      MTP = 'HO',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
      numCovar.1 = 1, tnum = 100000,
      R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2)
    pp_power <- pp1$D1indiv[2]
  }

  pp_power <- 0.97152

  vals <- test_sample_triad(pp_power, nbar = 50, J = 60, K = NULL,
                            seed = 22422422,
                            d_m = "d2.1_m2ff",
                            MTP = 'HO',
                            power.definition = 'D1indiv',
                            M = 3,
                            MDES = 0.125,
                            Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
                            numCovar.1 = 1, 
                            R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2)
  vals[1:2]
  warning_pattern(vals)

  expect_equal(60, vals$J, tol=0.1)
  expect_equal(50, vals$nbar, tol=0.1)
  expect_equal( warning_pattern(vals), c(FALSE,FALSE) )

})


test_that("testing of d2.1_m2fr two-tailed", {

  if ( FALSE ) {
    set.seed(8598)

    pp1 <- pump_power(
      d_m = "d2.1_m2fr",
      MTP = 'HO',
      nbar = 50,
      J = 60,
      M = 3,
      MDES = rep(0.125, 3),
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1, tnum = 100000,
      R2.1 = 0.1, ICC.2 = 0.05, omega.2 = 0.1, rho = 0.2)
    pp1
    pp_power <- pp1$D1indiv[2]
  }

  pp_power <- 0.92552

  vals <- test_sample_triad(pp_power, nbar = 50, J = 60, K = NULL,
                            seed = 22422422,
                            d_m = "d2.1_m2ff",
                            MTP = 'HO',
                            power.definition = 'D1indiv',
                            M = 3,
                            MDES = 0.125,
                            Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
                            numCovar.1 = 1,
                            R2.1 = 0.1, ICC.2 = 0.05, omega.2 = 0.1, rho = 0.2)
  vals[1:2]
  warning_pattern(vals)

  expect_equal(60, vals$J, tol = 0.1)
  expect_equal(50, vals$nbar, tol = 0.2)
  expect_equal( warning_pattern(vals), c(FALSE, FALSE) )

  mdes1 <-  pump_mdes(
    d_m = "d2.1_m2fr",
    MTP = 'HO',
    power.definition = 'D1indiv',
    target.power = pp_power,
    nbar = 50,
    J = 60,
    M = 3,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05, omega.2 = 0.1,
    rho = 0.2)
  expect_equal(mdes1$Adjusted.MDES, 0.125, tolerance = 0.1)

})


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -

# ------   d2.2_m2rc   -------

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # -

test_that("testing of d2.2_m2rc", {

  if ( FALSE ) {

    set.seed(8598)

    pp1 <- pump_power(
      d_m = "d2.2_m2rc",
      MTP = 'HO',
      nbar = 50,
      J = 20,
      M = 8,
      numZero = 5,
      MDES = 0.30,
      Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2,
      tnum = 100000)

    pp1
    pp_power <- pp1$min3[2]
    pp_power
  }

  set.seed( 4133333 )
  pp_power <- 0.66245

  vals <- test_sample_triad(pp_power, nbar = 50, J = 20, NULL, 2244323,
                            d_m = "d2.2_m2rc",
                            power.definition = "min3",
                            MTP = 'HO',
                            M = 8,
                            numZero = 5,
                            MDES = rep(0.30, 3),
                            Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
                            numCovar.1 = 1, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2 )
  
  vals
  expect_equal(20, vals$J, tol = 0.1 )
  expect_equal(50, vals$nbar, tol = 0.1 )
  expect_equal( warning_pattern(vals), c(FALSE, FALSE) )
  
  # cannot achieve target power with given parameters
  expect_warning(ss1 <- pump_sample(
      d_m = "d2.2_m2rc",
      MTP = 'BF',
      typesample = 'nbar',
      target.power = 0.8,
      power.definition = 'D1indiv',
      J = 20,
      M = 5,
      numZero = 0,
      MDES = 0.125,
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  ))
  
  expect_true(is.na(ss1$Sample.size))
  
  expect_warning(
    ss2 <- pump_sample(
      d_m = "d2.2_m2rc",
      MTP = 'HO',
      typesample = 'nbar',
      target.power = 0.8,
      power.definition = 'D1indiv',
      J = 20,
      M = 5,
      numZero = 0,
      MDES = 0.125,
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2 )
  )
 
  # even though we didn't converge, we get a value
  expect_true(!is.na(ss2$Sample.size))
  expect_true(ss2$Sample.size > 100000 )
  
  # can achieve target
  # checks power curve works for BF
  ss3 <- pump_sample(
      d_m = "d2.2_m2rc",
      MTP = 'BF',
      typesample = 'J',
      target.power = 0.8,
      power.definition = 'D1indiv',
      nbar = 200,
      M = 5,
      numZero = 0,
      MDES = 0.125,
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1, numCovar.2 = 1,
      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.3, rho = 0.2
  )
  
  expect_true(!is.null(plot_power_curve(ss3)))
  
})
