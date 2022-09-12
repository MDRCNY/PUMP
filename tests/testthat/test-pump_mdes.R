# library( PUMP )
# library( testthat )

test_that("pump_mdes runs for BF", {
  
  set.seed( 2424424 )

  pmdesB <- pump_mdes( d_m = "d2.1_m2fc",
                       MTP = "BF",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                       R2.1 = 0.1,ICC.2 = 0.05,
                       tnum = 1600,
                       rho = 0.4 )
  pmdesB
  expect_true( pmdesB$`Adjusted.MDES` > 0 )
  expect_true( abs(pmdesB$`D2indiv power` - 0.80) <  0.01 )
  
  # throw a warning correctly
  pmdesB <- expect_warning(pump_mdes( d_m = "d2.1_m2fc",
                       MTP = "BF",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                       R2.1 = 0.1,ICC.2 = 0.05,
                       tnum = 300,
                       rho = 0.4 ))
  pmdesB
  expect_true( pmdesB$`Adjusted.MDES` > 0 )
  expect_true( abs(pmdesB$`D2indiv power` - 0.80) <  0.01 )

  skip_on_cran()
  pmdesR <- expect_warning( pump_mdes( d_m = "d2.1_m2fc",
                       MTP = "None",
                       nbar = 200, J = 50,
                       power.definition = "D2indiv",
                       M = 3,
                       target.power = 0.80, tol = 0.01,
                       Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                       R2.1 = 0.1,ICC.2 = 0.05,
                       tnum = 300,
                       rho = 0.4 ) )
  expect_true( pmdesR$`Adjusted.MDES` > 0 )
  expect_true( abs(pmdesR$`D2indiv power` - 0.80) <  0.01 )
  expect_true( pmdesR$`Adjusted.MDES` < pmdesB$`Adjusted.MDES`)

  set.seed( 14444444 )
  pmdesBmin <- pump_mdes(
                      d_m = "d2.1_m2fc",
                      MTP = "BF",
                      nbar = 200, J = 50,
                      power.definition = "min1",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                      R2.1 = 0.1, ICC.2 = 0.05,
                      tnum = 1600,
                      rho = 0.4 )
  pmdesBmin
  expect_true( pmdesBmin$`Adjusted.MDES` < pmdesR$`Adjusted.MDES` )
  expect_true( abs( pmdesBmin$`min1 power` - 0.80) <  0.01 )

  set.seed( 444224 )
  pmdes_comp <- pump_mdes( d_m = "d2.1_m2fc",
                           MTP = "BF",
                           nbar = 200, J = 50,
                           power.definition = "complete",
                           M = 3,
                           target.power = 0.80, tol = 0.02,
                           Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                           R2.1 = 0.1, ICC.2 = 0.05,
                           tnum = 400,
                           rho = 0.4 )
  pmdes_comp
  expect_true( pmdes_comp$`Adjusted.MDES` > pmdesB$`Adjusted.MDES` )
  expect_true( abs( pmdes_comp$`complete power` - 0.80) <  0.02 )

  sp <- search_path( pmdes_comp )
  expect_true( is.data.frame(sp) )
  expect_true( max( sp$w ) == 400*4 )
  
  ppBcomp <- pump_power(
                   d_m = "d2.1_m2fc",
                   MTP = "BF",
                   MDES = rep( pmdes_comp$`Adjusted.MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                   R2.1 = 0.1, ICC.2 = 0.05,
                   rho = 0.4 )

  ppBcomp
  expect_true( abs( ppBcomp[2,"complete"] - 0.80 ) <= 0.02)

  pmdesBmin$mdes.results
  ppBmin <- pump_power(
                   d_m = "d2.1_m2fc",
                   MTP = "BF",
                   MDES = rep( pmdesBmin$`Adjusted.MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                   R2.1 = 0.1,ICC.2 = 0.05,
                   rho = 0.4 )
  ppBmin
  expect_true( abs( ppBmin[2,"min1"] - 0.80 ) <= 0.02)

})


test_that("pump_mdes runs for D1indiv, HO", {
    
  skip_on_cran()

  set.seed( 1010101 )
  pmdes <- pump_mdes( d_m = "d2.1_m2fc",
                      MTP = "HO",
                      nbar = 200, J = 50,
                      power.definition = "D1indiv",
                      M = 3,
                      target.power = 0.80, tol = 0.01,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                      R2.1 = 0.1, ICC.2 = 0.05,
                      tnum = 1600,
                      rho = 0.4 )

  pmdes
  expect_true( pmdes$`Adjusted.MDES` > 0 )
  expect_true( abs( pmdes$`D1indiv power` - 0.80) <  0.01 )

  pp <- pump_power( d_m = "d2.1_m2fc",
                   MTP = "HO",
                   MDES = rep( pmdes$`Adjusted.MDES`, 3 ),
                   nbar = 200, J = 50,
                   M = 3,
                   Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                   R2.1 = 0.1, ICC.2 = 0.05,
                   rho = 0.4 )
  pp
  expect_true( abs( pp[2,2] - 0.80 ) <= 0.02)
})


test_that("pump_mdes runs for d1.1_m1c", {
    
  skip_on_cran()

  set.seed( 10130103 )
  R2.1 <- 0.61
  pmdes <- pump_mdes(d_m = "d1.1_m1c", MTP = "HO",
                     target.power = 0.80, power.definition = "min1", tol = 0.01,
                     R2.1 = R2.1, numCovar.1 = 1, J = 1,
                     tnum = 1600,
                     M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5)
  pmdes

  ES <- pmdes$`Adjusted.MDES`
  ppow <- pump_power(d_m = "d1.1_m1c", MTP = "HO", MDES = ES,
             R2.1 = R2.1, numCovar.1 = 1,
             M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5 )
  ppow

  expect_true( !is.null( pmdes ) )
  expect_true( abs( ppow$min1[[2]]- 0.80 ) <= 0.02 )
} )


test_that("No adjustment", {
    
    skip_on_cran()

    pmdes <- expect_warning( pump_mdes( d_m = "d2.1_m2fc",
                        MTP = "None",
                        nbar = 200, J = 50,
                        power.definition = "D1indiv",
                        M = 3,
                        target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        tnum = 1600,
                        rho = 0.4 ) )

    expect_error(expect_warning(pmdes <- pump_mdes( d_m = "d2.1_m2fc",
                        MTP = "None",
                        nbar = 200, J = 50,
                        power.definition = "min2",
                        M = 3,
                        target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                        R2.1 = 0.1,ICC.2 = 0.05,
                        tnum = 1600,
                        rho = 0.4 )))


    pmdes <- pump_mdes( d_m = "d2.1_m2fc",
                        MTP = "HO",
                        nbar = 200, J = 50,
                        power.definition = "D1indiv",
                        M = 3,
                        target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        tnum = 1600,
                        rho = 0.4 )
})

test_that("power definitions", {
    
  skip_on_cran()

  pmdes <- pump_mdes( d_m = "d2.1_m2fc",
                      MTP = "HO",
                      nbar = 200, J = 50,
                      power.definition = "indiv.mean",
                      M = 3,
                      target.power = 0.80,
                      Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
                      R2.1 = 0.1, ICC.2 = 0.05,
                      tnum = 400,
                      rho = 0.4 )
  expect_true(!is.null(pmdes))
})


test_that( "errors out when providing MDES", {
    
  skip_on_cran()
    
  expect_error(pmdes <- pump_mdes(
    d_m = "d2.1_m2fc",
    MDES = rep(0.2, 5),
    MTP = "HO",
    nbar = 200, J = 50,
    power.definition = "indiv.mean",
    M = 3,
    target.power = 0.80, tol = 0.01,
    Tbar = 0.50, alpha = 0.05, numCovar.1 = 5,
    R2.1 = 0.1, ICC.2 = 0.05,
    tnum = 300,
    rho = 0.4 ))
  

})

test_that( "different values for different outcomes", {
    
  skip_on_cran()

  set.seed(03443)

  pow <- pump_power(
    d_m = "d2.1_m2fc",
    MTP = "HO",
    J = 20,
    nbar = 200,
    M = 3,
    MDES = 0.05,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5,
    R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
    rho = 0.4 )

  # sanity check: higher ICC means higher power
  expect_true(pow$D2indiv[1] > pow$D1indiv[1])
  expect_true(pow$D3indiv[1] > pow$D2indiv[1])

  mdes1 <- pump_mdes(
    d_m = "d2.1_m2fc",
    MTP = "HO",
    target.power = 0.8,
    power.definition = 'D1indiv',
    J = 20,
    nbar = 200,
    M = 3,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5,
    R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
    rho = 0.4
  )

  mdes2 <- pump_mdes(
    d_m = "d2.1_m2fc",
    MTP = "HO",
    target.power = 0.8,
    power.definition = 'D2indiv',
    J = 20,
    nbar = 200,
    M = 3,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5,
    R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
    rho = 0.4
  )

  mdes3 <- pump_mdes(
    d_m = "d2.1_m2fc",
    MTP = "HO",
    target.power = 0.8,
    power.definition = 'D3indiv',
    J = 20,
    nbar = 200,
    M = 3,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5,
    R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
    rho = 0.4
  )

  # for same target power, we should have a smaller MDES for larger ICC
  expect_true(mdes1$Adjusted.MDES > mdes2$Adjusted.MDES)
  expect_true(mdes2$Adjusted.MDES > mdes3$Adjusted.MDES)

})


test_that("M > 1 with MTP None", {

    skip_on_cran()
    
    pmdes <- expect_warning(pump_mdes(
                          d_m = "d2.1_m2fc",
                          target.power = 0.8,
                          power.definition = 'D1indiv',
                          MTP = "None",
                          M = 3,
                          J = 3, # number of schools/block
                          nbar = 258,
                          Tbar = 0.50, # prop Tx
                          alpha = 0.05, # significance level
                          numCovar.1 = 5,
                          R2.1 = 0.1,
                          ICC.2 = 0.05,
                          rho = 0.4
    ))
    expect_true( nrow( pmdes ) == 1 )
    
    expect_error(expect_warning(pump_mdes( d_m = "d2.1_m2fc",
                         target.power = 0.8,
                         power.definition = "complete",
                         MTP = "None",
                         M = 3,
                         J = 3, # number of schools/block
                         nbar = 258,
                         Tbar = 0.50, # prop Tx
                         alpha = 0.05, # significance level
                         numCovar.1 = 5,
                         R2.1 = 0.1,
                         ICC.2 = 0.05,
                         rho = 0.4
    )))
})
