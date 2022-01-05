# library( PUMP )
# library( testthat )
library(PowerUpR)

default.tol <- 0.01

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d2.1_m2fc -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fc", {

  sink("sink.txt")

  powerup.power <- PowerUpR::power.bira2c1(
    es = 0.125,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g1 = 1,
    r21 = 0.1,
    n = 50,
    J = 30
  )

  pump.power <- pump_power(
    d_m = "d2.1_m2fc",
    MTP = 'None',
    nbar = 50,
    J = 30,
    M = 1,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0, rho = 0,
    tnum = 100000
  )

  expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)
  
  powerup.power <- PowerUpR::power.bira2c1(
      es = 0.125,
      alpha = 0.05,
      two.tailed = TRUE,
      p = 0.5,
      g1 = 1,
      r21 = 0.1,
      n = 50,
      J = 30
  )
  
  pump.power <- pump_power(
      d_m = "d2.1_m2fc",
      MTP = 'None',
      nbar = 50,
      J = 30,
      M = 1,
      MDES = 0.125,
      Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
      numCovar.1 = 1,
      R2.1 = 0.1, ICC.2 = 0, rho = 0,
      tnum = 100000
  )
  
  expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)

  sink()
  file.remove("sink.txt")
  
})


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d2.1_m2fr -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d2.1_m2fr one-tailed", {

  sink("sink.txt")

  powerup.power <- expect_warning(PowerUpR::power.bira2r1(
    es = 0.125,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    omega2 = 0.1,
    r21 = 0.3,
    r2t2 = 0,
    n = 50,
    J = 30
  ))

  pump.power <- pump_power(
    d_m = "d2.1_m2fr",
    MTP = 'None',
    nbar = 50,
    J = 30,
    M = 1,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
    omega.2 = 0.1,
    tnum = 100000
  )

  expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)

  powerup.mdes <- expect_warning(PowerUpR::mdes.bira2r1(
    power = 0.8,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    omega2 = 0.1,
    r21 = 0.3,
    r2t2 = 0,
    n = 50,
    J = 30
  ))

  pump.mdes <- pump_mdes(
    target.power = 0.8,
    power.definition = 'D1indiv',
    d_m = "d2.1_m2fr",
    MTP = 'None',
    nbar = 50,
    J = 30,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
    omega.2 = 0.1
  )

  expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)


  powerup.ss <- expect_warning(PowerUpR::mrss.bira2r1(
    power = 0.8,
    es = 0.125,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    omega2 = 0.1,
    r21 = 0.3,
    r2t2 = 0,
    n = 50
  ))

  pump.ss <- pump_sample(
    target.power = 0.8,
    power.definition = 'D1indiv',
    typesample = 'J',
    d_m = "d2.1_m2fr",
    MTP = 'None',
    MDES = 0.125,
    nbar = 50,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
    omega.2 = 0.1
  )

  expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)

  sink()
  file.remove("sink.txt")

})


test_that("testing of d2.1_m2fr two-tailed", {

  sink("sink.txt")

  powerup.power <- expect_warning(PowerUpR::power.bira2r1(
    es = 0.125,
    alpha = 0.05,
    two.tailed = TRUE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    omega2 = 0.1,
    r21 = 0.3,
    r2t2 = 0,
    n = 50,
    J = 30
  ))

  pump.power <- pump_power(
    d_m = "d2.1_m2fr",
    MTP = 'None',
    nbar = 50,
    J = 30,
    M = 1,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
    omega.2 = 0.1,
    tnum = 100000
  )

  expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)

  powerup.mdes <- expect_warning(PowerUpR::mdes.bira2r1(
    power = 0.8,
    alpha = 0.05,
    two.tailed = TRUE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    omega2 = 0.1,
    r21 = 0.3,
    r2t2 = 0,
    n = 50,
    J = 30
  ))

  pump.mdes <- pump_mdes(
    target.power = 0.8,
    power.definition = 'D1indiv',
    d_m = "d2.1_m2fr",
    MTP = 'None',
    nbar = 50,
    J = 30,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
    omega.2 = 0.1
  )

  expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)


  powerup.ss <- expect_warning(PowerUpR::mrss.bira2r1(
    power = 0.8,
    es = 0.125,
    alpha = 0.05,
    two.tailed = TRUE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    omega2 = 0.1,
    r21 = 0.3,
    r2t2 = 0,
    n = 50
  ))

  pump.ss <- pump_sample(
    target.power = 0.8,
    power.definition = 'D1indiv',
    typesample = 'J',
    d_m = "d2.1_m2fr",
    MTP = 'None',
    MDES = 0.125,
    nbar = 50,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.3, ICC.2 = 0.05, rho = 0,
    omega.2 = 0.1
  )

  expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)

  sink()
  file.remove("sink.txt")

})


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ------------- d3.2_m3ff2rc -------------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

test_that("testing of d3.2_m3ff2rc one-tailed", {

  sink("sink.txt")

  powerup.power <- PowerUpR::power.bcra3f2(
    es = 0.125,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    r21 = 0.3,
    r22 = 0.3,
    n = 50,
    J = 30,
    K = 10
  )

  pump.power <- pump_power(
    d_m = "d3.2_m3ff2rc",
    MTP = 'None',
    nbar = 50,
    J = 30,
    K = 10,
    M = 1,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
    tnum = 100000
  )

  expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)

  powerup.mdes <- PowerUpR::mdes.bcra3f2(
    power = 0.8,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    r21 = 0.3,
    r22 = 0.3,
    n = 50,
    J = 30,
    K = 10
  )

  pump.mdes <- pump_mdes(
    target.power = 0.8,
    power.definition = 'D1indiv',
    d_m = "d3.2_m3ff2rc",
    MTP = 'None',
    nbar = 50,
    J = 30,
    K = 10,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
  )

  expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)


  powerup.ss <- PowerUpR::mrss.bcra3f2(
    power = 0.8,
    es = 0.125,
    alpha = 0.05,
    two.tailed = FALSE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    r21 = 0.3,
    r22 = 0.3,
    n = 50,
    J = 30,
    K = 10
  )

  pump.ss <- pump_sample(
    target.power = 0.8,
    power.definition = 'D1indiv',
    typesample = 'J',
    d_m = "d3.2_m3ff2rc",
    MTP = 'None',
    MDES = 0.125,
    nbar = 50,
    K = 10,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
    tol = 0.005
  )

  print(pump.ss$Sample.size)
  print(powerup.ss$J)
  expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)

  sink()
  file.remove("sink.txt")

})

test_that("testing of d3.2_m3ff2rc two-tailed", {

  sink("sink.txt")

  powerup.power <- PowerUpR::power.bcra3f2(
    es = 0.125,
    alpha = 0.05,
    two.tailed = TRUE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    r21 = 0.3,
    r22 = 0.3,
    n = 50,
    J = 30,
    K = 10
  )

  pump.power <- pump_power(
    d_m = "d3.2_m3ff2rc",
    MTP = 'None',
    nbar = 50,
    J = 30,
    K = 10,
    M = 1,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
    tnum = 100000
  )

  expect_equal(pump.power$D1indiv[1], powerup.power$power, tol = default.tol)

  powerup.mdes <- PowerUpR::mdes.bcra3f2(
    power = 0.8,
    alpha = 0.05,
    two.tailed = TRUE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    r21 = 0.3,
    r22 = 0.3,
    n = 50,
    J = 30,
    K = 10
  )

  pump.mdes <- pump_mdes(
    target.power = 0.8,
    power.definition = 'D1indiv',
    d_m = "d3.2_m3ff2rc",
    MTP = 'None',
    nbar = 50,
    J = 30,
    K = 10,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
  )

  expect_equal(pump.mdes$Adjusted.MDES, powerup.mdes$mdes[1], tol = default.tol)


  powerup.ss <- PowerUpR::mrss.bcra3f2(
    power = 0.8,
    es = 0.125,
    alpha = 0.05,
    two.tailed = TRUE,
    p = 0.5,
    g2 = 1,
    rho2 = 0.05,
    r21 = 0.3,
    r22 = 0.3,
    n = 50,
    J = 30,
    K = 10
  )

  pump.ss <- pump_sample(
    target.power = 0.8,
    power.definition = 'D1indiv',
    typesample = 'J',
    d_m = "d3.2_m3ff2rc",
    MTP = 'None',
    MDES = 0.125,
    nbar = 50,
    K = 10,
    M = 1,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.3, R2.2 = 0.3, ICC.2 = 0.05, rho = 0,
    tol = 0.005
  )
  
  expect_equal(pump.ss$Sample.size, powerup.ss$J, tol = default.tol)

  sink()
  file.remove("sink.txt")

})

