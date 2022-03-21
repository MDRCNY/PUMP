# library( PUMP )
# library( testthat )
default.tnum <- 200

test_that("pump_power works without crashing", {

  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 21, # number RA blocks
                    nbar = 258,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, # how correlated outcomes are
                    tnum = 200
  )
  pp
  expect_true( is.pumpresult( pp ) )

  expect_equal( dim( pp ), c(2,8) )

  expect_true( is.na( pp[1,6] ) )
  expect_true( is.na( pp[1,7] ) )
  expect_true( is.na( pp[1,8] ) )
})


test_that( "pump_power handles small MDES right", {
  
  skip_on_cran()
    
  set.seed(58554343)
  
  aa <- pump_power(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    MDES = 0.000001,
    J = 60,
    nbar = 50,
    M = 3,
    Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05,
    rho = 0.2,
    tnum = default.tnum)
  aa
  
  expect_equal( 0.05, aa$indiv.mean[[1]], tolerance = 0.01 )
  
  aa <- pump_power(
    d_m = "d2.1_m2fc",
    MTP = 'HO',
    MDES = 0.000001,
    J = 60,
    nbar = 50,
    M = 3,
    Tbar = 0.5, alpha = 0.05, two.tailed = FALSE,
    numCovar.1 = 1,
    R2.1 = 0.1, ICC.2 = 0.05,
    rho = 0.2,
    tnum = default.tnum)
  aa
  expect_equal( 0.05, aa$indiv.mean[[1]], tolerance = 0.01 )
  
})



test_that("pump_power long.table", {

  skip_on_cran()
    
  set.seed(9515)

  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = "BF",
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 3, # number of schools/block
                    K = 21, # number RA blocks
                    nbar = 258,
                    Tbar = 0.50, # prop Tx
                    alpha = 0.05, # significance level
                    numCovar.1 = 5, numCovar.2 = 3,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.4,
                    rho = 0.4, # how correlated outcomes are
                    tnum = 200,
                    long.table = TRUE
  )
  pp
  expect_true( is.pumpresult( pp ) )
  # expect_true( is.numeric( pp$BF ) )

  expect_equal( dim( pp ), c(7,3) )

  expect_true( is.na( pp$None[5] ) )
  expect_true( is.na( pp$None[6] ) )
  expect_true( is.na( pp$None[7] ) )
  expect_true( pp$BF[4] < pp$BF[5] )
  expect_true( pp$BF[5] > pp$BF[6] )
  expect_true( pp$BF[6] > pp$BF[7] )
  expect_true( pp$BF[4] > pp$BF[7] )
  expect_true( all ( ( pp$None >= pp$BF )[1:4] ) )
})





test_that("skipping level three inputs for level 2 works", {
    
  skip_on_cran()

  expect_warning(pp <- pump_power(   d_m = "d2.1_m2fc",
                                     MTP = "BF",
                                     MDES = rep( 0.10, 3 ),
                                     M = 3,
                                     J = 3, # number of schools/block
                                     K = 21, # number RA blocks
                                     nbar = 258,
                                     Tbar = 0.50, # prop Tx
                                     alpha = 0.05, # significance level
                                     numCovar.1 = 5, numCovar.2 = 3,
                                     R2.1 = 0.1, R2.2 = 0.7,
                                     ICC.2 = 0.05,
                                     rho = 0.4, tnum = 200
  ))

  expect_equal( dim( pp ), c(2,8) )

  expect_true( pp[2,"min1"] >= pp[2,"D1indiv"] )
})

test_that("having no covariates is fine", {

  skip_on_cran()
    
  pp <- pump_power(   d_m = "d2.1_m2fc",
                      MTP = "BF",
                      MDES = rep( 0.10, 3 ),
                      M = 3,
                      J = 3, # number of schools/block
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      ICC.2 = 0.05,
                      rho = 0.4, tnum = 200
  )

  expect_equal( dim( pp ), c(2,8) )
})




test_that("pump_power works with multiple MTP", {
  pp <- pump_power( d_m = "d3.2_m3ff2rc",
                    MTP = c( "HO", "BF" ),
                    MDES = rep( 0.10, 3 ),
                    M = 3,
                    J = 20, # number of schools/block
                    K = 21, # number RA blocks
                    nbar = 12,
                    Tbar = 0.20, # prop Tx
                    alpha = 0.15, # significance level
                    numCovar.1 = 4, numCovar.2 = 1,
                    R2.1 = 0.1, R2.2 = 0.7,
                    ICC.2 = 0.05, ICC.3 = 0.1,
                    rho = 0.4, tnum = 200
  ) # how correlated outcomes are
  expect_equal( dim( pp ), c(3,8) )

})



test_that("M = 1 runs successfully", {

  skip_on_cran()
    
  pp <- pump_power(   d_m = "d2.1_m2fc",
                      MTP = "None",
                      MDES = 0.10,
                      M = 1,
                      J = 3, # number of schools/block
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5,
                      R2.1 = 0.1,
                      ICC.2 = 0.05,
                      rho = 0.4, tnum = 200
  )
  expect_true( nrow( pp ) == 1 )
  
  pp <- pump_power(   d_m = "d2.1_m2fc",
                      MTP = "None",
                      MDES = 0.10,
                      M = 1,
                      J = 3, # number of schools/block
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5,
                      R2.1 = 0.1,
                      ICC.2 = 0.05,
                      rho = 0.4, tnum = 200
  )
  expect_true( nrow( pp ) == 1 )
  expect_true( ncol( pp ) == 2 )
})

test_that("J = 1 runs successfully", {

    skip_on_cran()

    expect_warning(pp <- pump_power(   d_m = "d2.1_m2fc",
                        MTP = "None",
                        MDES = 0.10,
                        M = 1,
                        J = 1, # number of schools/block
                        nbar = 258,
                        Tbar = 0.50, # prop Tx
                        alpha = 0.05, # significance level
                        numCovar.1 = 5, numCovar.2 = 3,
                        R2.1 = 0.1, R2.2 = 0.7,
                        ICC.2 = 0.05,
                        rho = 0.4, tnum = 200
    ))
    expect_true( nrow( pp ) == 1 )
})

test_that("K = 1 runs successfully", {

    skip_on_cran()
    expect_warning(pp <- pump_power( d_m = "d3.2_m3ff2rc",
                      MTP = "BF",
                      MDES = rep( 0.10, 3 ),
                      M = 3,
                      J = 10, # number of schools/block
                      K = 1, # number RA blocks
                      nbar = 150,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 0, numCovar.2 = 0,
                      R2.1 = 0, R2.2 = 0,
                      ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0.4, tnum = 200
    ))
    pp
})

test_that("unblocked d_ms", {

  skip_on_cran()

  pp <- pump_power(   d_m = "d1.1_m1c",
                      MTP = "BF",
                      MDES = rep( 0.50, 3 ),
                      M = 3,
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5,
                      R2.1 = 0.1,
                      rho = 0.4, tnum = 200
  )


  ES <- log( 2 ) / 0.66
  ES
  R2.2 <- 0.6102
  pump_power(d_m = "d1.1_m1c", MTP = "HO", MDES = ES,
             M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5 )


  expect_true( nrow( pp ) == 2 )

})

test_that("Correct MTP parameter validation.", {

    skip_on_cran()

    pp <- expect_warning(pump_power(   d_m = "d2.1_m2fc",
                                       MTP = "None",
                                       MDES = rep( 0.10, 3 ),
                                       M = 3,
                                       J = 3, # number of schools/block
                                       nbar = 258,
                                       Tbar = 0.50, # prop Tx
                                       alpha = 0.05, # significance level
                                       numCovar.1 = 5,
                                       R2.1 = 0.1,
                                       ICC.2 = 0.05,
                                       rho = 0.4, tnum = 200
    ))

    # no MTP provided
    expect_error(pp <- pump_power( d_m = "d2.1_m2fc",
                                     MDES = 0.1,
                                     M = 3,
                                     J = 3, # number of schools/block
                                     nbar = 258,
                                     Tbar = 0.50, # prop Tx
                                     alpha = 0.05, # significance level
                                     numCovar.1 = 5,
                                     R2.1 = 0.1,
                                     ICC.2 = 0.05,
                                     rho = 0.4, tnum = 200
    ))

    # no MTP provided with single outcome is fine.
    pp <- pump_power( d_m = "d2.1_m2fc",
                                   MDES = 0.1,
                                   M = 1,
                                   J = 3, # number of schools/block
                                   nbar = 258,
                                   Tbar = 0.50, # prop Tx
                                   alpha = 0.05, # significance level
                                   numCovar.1 = 5,
                                   R2.1 = 0.1,
                                   ICC.2 = 0.05,
                                   rho = 0.4, tnum = 200
    )
    expect_true( nrow( pp ) == 1 )
})


test_that("numZero has expected behavior", {
    
  skip_on_cran()

  pp <- pump_power( d_m = "d2.2_m2rc",
                    MTP = "BF",
                    J = 10,
                    M = 5,
                    nbar = 100,
                    MDES = rep( 0.2, 2 ),
                    ICC.2 = 0.1, R2.2 = 0.3, numCovar.2 = 1,
                    numZero = 3,
                    Tbar = 0.50,
                    rho = 1,
                    tnum = 200)



  expect_error(pp <- pump_power( d_m = "d2.2_m2rc",
                    MTP = "BF",
                    J = 10,
                    M = 4,
                    nbar = 100,
                    MDES = rep( 0.2, 2 ),
                    ICC.2 = 0.1, R2.2 = 0.3, numCovar.2 = 1,
                    numZero = 3,
                    Tbar = 0.50,
                    rho = 1,
                    tnum = 200)
  )

})

test_that("do not report invalid power values", {

  skip_on_cran()

  pp <- pump_power( d_m = "d2.2_m2rc",
                    MTP = "BF",
                    J = 10,
                    M = 3,
                    nbar = 100,
                    MDES = 0.2,
                    Tbar = 0.50,
                    ICC.2 = 0.1,
                    R2.2 = 0.1, numCovar.2 = 1,
                    rho = 0.4,
                    tnum = 200)

  expect_true(is.na(pp$min1[1]))
  expect_true(is.na(pp$min2[1]))
  expect_true(is.na(pp$complete[1]))
})

test_that("M > 1 with MTP None runs successfully", {

    skip_on_cran()
    
    pp <- expect_warning(pump_power(   d_m = "d2.1_m2fc",
                        MTP = "None",
                        MDES = 0.10,
                        M = 3,
                        J = 3, # number of schools/block
                        nbar = 258,
                        Tbar = 0.50, # prop Tx
                        alpha = 0.05, # significance level
                        numCovar.1 = 5,
                        R2.1 = 0.1,
                        ICC.2 = 0.05,
                        rho = 0.4, tnum = 200
    ))
    expect_true( nrow( pp ) == 1 )
    
    pp <- pump_power(   d_m = "d2.1_m2fc",
                                       MTP = c("BF", "None"),
                                       MDES = 0.10,
                                       M = 3,
                                       J = 3, # number of schools/block
                                       nbar = 258,
                                       Tbar = 0.50, # prop Tx
                                       alpha = 0.05, # significance level
                                       numCovar.1 = 5,
                                       R2.1 = 0.1,
                                       ICC.2 = 0.05,
                                       rho = 0.4, tnum = 200
    )
    expect_true( nrow( pp ) == 2)
})

test_that("zero MDES values", {
    
    skip_on_cran()
    
    # zero in middle of vector
    pp <- pump_power(   d_m = "d2.1_m2fc",
                                       MTP = "HO",
                                       MDES = c(0.1, 0, 0.1),
                                       M = 3,
                                       J = 3, # number of schools/block
                                       nbar = 258,
                                       Tbar = 0.50, # prop Tx
                                       alpha = 0.05, # significance level
                                       numCovar.1 = 5,
                                       R2.1 = 0.1,
                                       ICC.2 = 0.05,
                                       rho = 0.4, tnum = 200
    )
    
    expect_true(
        all(colnames(pp) ==
        c('MTP', 'D1indiv', 'D3indiv', 'indiv.mean', 'min1', 'min2', 'complete'))
    )
    
    # all zero, don't drop zero outcomes
    pp <- pump_power(   d_m = "d2.1_m2fc",
                        MTP = "HO",
                        MDES = 0,
                        M = 3,
                        J = 3, # number of schools/block
                        nbar = 258,
                        Tbar = 0.50, # prop Tx
                        alpha = 0.05, # significance level
                        numCovar.1 = 5,
                        R2.1 = 0.1,
                        ICC.2 = 0.05,
                        rho = 0.4, tnum = 200,
                        drop.zero.outcomes = FALSE
    )
    expect_true(nrow(pp) == 2)
    
    
    # all zero, do drop zero outcomes
    pp <- pump_power(   d_m = "d2.1_m2fc",
                        MTP = "HO",
                        MDES = 0,
                        M = 3,
                        J = 3, # number of schools/block
                        nbar = 258,
                        Tbar = 0.50, # prop Tx
                        alpha = 0.05, # significance level
                        numCovar.1 = 5,
                        R2.1 = 0.1,
                        ICC.2 = 0.05,
                        rho = 0.4, tnum = 200,
                        drop.zero.outcomes = TRUE
    )
    
    expect_true(is.na(pp$D1indiv[1]) && is.na(pp$D1indiv[2]))
})

test_that("different correlations", {
    
    skip_on_cran()
    
    pp.rhomin <- pump_power( d_m = "d2.2_m2rc",
                             MTP = "BF",
                             J = 10,
                             M = 4,
                             nbar = 100,
                             MDES = rep( 0.2, 4 ),
                             Tbar = 0.50, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.5, ICC.2 = 0.05,
                             rho = 0, tnum = 1000)
    
    pp.rhomed <- pump_power(   d_m = "d2.2_m2rc",
                               MTP = "BF",
                               J = 10,
                               M = 4,
                               nbar = 100,
                               MDES = rep( 0.2, 4 ),
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.5, ICC.2 = 0.05,
                               rho = 0.4, tnum = 1000)
    
    pp.rhomax <- pump_power(   d_m = "d2.2_m2rc",
                               MTP = "BF",
                               J = 10,
                               M = 4,
                               nbar = 100,
                               MDES = rep( 0.2, 4 ),
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.5, ICC.2 = 0.05,
                               rho = 1, tnum = 1000)
    
    pp.rhomin
    pp.rhomed
    pp.rhomax
    # lower correlation means higher min1 power (more chances to hit one out of the park)
    expect_true( pp.rhomin$min1[2] > pp.rhomed$min1[2] )
    expect_true( pp.rhomed$min1[2] > pp.rhomax$min1[2]  )
    
    # complete power is the reverse
    expect_true( pp.rhomin$complete[2] < pp.rhomed$complete[2] )
    expect_true( pp.rhomed$complete[2] < pp.rhomax$complete[2] )
    
})

test_that("testing against Porter 2017", {
    
    skip_on_cran()
    
    # NOTE: compare to Figure 1, setting with:
    # all outcomes have effects
    # 3 outcomes
    
    set.seed(13434)
    pp1 <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = c('BF', 'HO', 'BH', 'WY-SS', 'WY-SD'),
        nbar = 50,
        J = 20,
        M = 3,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, tnum = 1000, B = 1000,
        R2.1 = 0.5, ICC.2 = 0, rho = 0.2)
    
    expect_equal(0.65, pp1$D1indiv[2], tol = 0.01)
    expect_equal(0.71, pp1$D1indiv[3], tol = 0.01)
    expect_equal(0.76, pp1$D1indiv[4], tol = 0.01)
    expect_equal(0.67, pp1$D1indiv[5], tol = 0.01)
    expect_equal(0.75, pp1$D1indiv[6], tol = 0.01)
    
    set.seed(13434)
    pp2 <- pump_power(
        d_m = "d2.1_m2fc",
        MTP = c('BF', 'HO', 'BH', 'WY-SS', 'WY-SD'),
        nbar = 50,
        J = 20,
        M = 3,
        MDES = 0.125,
        Tbar = 0.5, alpha = 0.05, two.tailed = TRUE,
        numCovar.1 = 1, tnum = 1000, B = 1000,
        R2.1 = 0.5, ICC.2 = 0, rho = 0.8)
    
    expect_equal(0.65, pp2$D1indiv[2], tol = 0.01)
    expect_equal(0.71, pp2$D1indiv[3], tol = 0.01)
    expect_equal(0.75, pp2$D1indiv[4], tol = 0.01)
    expect_equal(0.71, pp2$D1indiv[5], tol = 0.01)
    expect_equal(0.76, pp2$D1indiv[6], tol = 0.01)
    
})


test_that("MTP behavior", {
    
    skip_on_cran()
    
    set.seed(13434)
    pp <- pump_power( d_m = "d3.2_m3ff2rc",
                      MTP = c("BF", "HO", "BH", "WY-SS", "WY-SD"),
                      MDES = c(0.025, 0.05, 0.1, 0.15, 0.2),
                      M = 5,
                      J = 3, # number of schools/block
                      K = 10, # number RA blocks
                      nbar = 258,
                      Tbar = 0.50, # prop Tx
                      alpha = 0.05, # significance level
                      numCovar.1 = 5, numCovar.2 = 3,
                      R2.1 = 0.1, R2.2 = 0.7,
                      ICC.2 = 0.05, ICC.3 = 0.4,
                      rho = 0.4, # how correlated outcomes are
                      tnum = 1000, B = 1000
    )
    pp
    
    # for biggest effect, BF same as HO
    expect_equal(pp$D5indiv[pp$MTP == 'BF'], pp$D5indiv[pp$MTP == 'HO'], tol = 0.1)
    # for smaller effects, HO more powerful than Bonf
    expect_true(pp$D1indiv[pp$MTP == 'BF'] <pp$D1indiv[pp$MTP == 'HO'])
    
    # BH least conservative (for small effects)
    expect_true(pp$D1indiv[pp$MTP == 'BF'] < pp$D1indiv[pp$MTP == 'BH'])
    
    # WY more powerful than BF
    expect_true(pp$D5indiv[pp$MTP == 'BF'] < pp$D5indiv[pp$MTP == 'WY-SS'])
    # for biggest effect, WY-SD similar to WY-SS for largest outcome
    expect_equal(pp$D1indiv[pp$MTP == 'WY-SS'], pp$D1indiv[pp$MTP == 'WY-SD'], tol = 0.1)
    # for smaller effects, WY-SD more powerful than WY-SS for largest outcome
    expect_true(pp$D1indiv[pp$MTP == 'WY-SS'] < pp$D1indiv[pp$MTP == 'WY-SD'])
    
})

test_that("parallel WY", {
    
    skip_on_cran()
    
    pp <- expect_warning(pump_power( d_m = "d3.2_m3ff2rc",
                                     MTP = c("WY-SD"),
                                     MDES = c(0.025, 0.05, 0.1, 0.15, 0.2),
                                     M = 5,
                                     J = 3, # number of schools/block
                                     K = 10, # number RA blocks
                                     nbar = 258,
                                     Tbar = 0.50, # prop Tx
                                     alpha = 0.05, # significance level
                                     numCovar.1 = 5, numCovar.2 = 3,
                                     R2.1 = 0.1, R2.2 = 0.7,
                                     ICC.2 = 0.05, ICC.3 = 0.4,
                                     rho = 0.4, # how correlated outcomes are
                                     tnum = 200, B = 100,
                                     parallel.WY.cores = 1,
                                     verbose = TRUE
    ))
    expect_true(nrow(pp) == 2)
    
    
    pp <- expect_warning(pump_power( d_m = "d3.2_m3ff2rc",
                                     MTP = c("WY-SD"),
                                     MDES = c(0.025, 0.05, 0.1, 0.15, 0.2),
                                     M = 5,
                                     J = 3, # number of schools/block
                                     K = 10, # number RA blocks
                                     nbar = 258,
                                     Tbar = 0.50, # prop Tx
                                     alpha = 0.05, # significance level
                                     numCovar.1 = 5, numCovar.2 = 3,
                                     R2.1 = 0.1, R2.2 = 0.7,
                                     ICC.2 = 0.05, ICC.3 = 0.4,
                                     rho = 0.4, # how correlated outcomes are
                                     tnum = 200, B = 100,
                                     parallel.WY.cores = 2
    ))
    expect_true(nrow(pp) == 2)
    
})


