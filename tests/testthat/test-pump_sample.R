# library( PUMP )
# library( testthat )

default.tnum <- 1000

test_that("calc_nbar works", {

  nbar <- PUMP:::calc_nbar(  d_m = "d2.2_m2rc",
                            MT = 2.8,
                            MDES = 0.20,
                            J = 5,
                            Tbar = 0.25,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  nbar
  expect_true( is.na( nbar ) )

  nbar <- calc_nbar(  d_m = "d2.2_m2rc",
                            MT = 2.8,
                            MDES = 0.20,
                            J = 305,
                            Tbar = 0.25,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  nbar
  expect_true( !is.na( nbar ) )

} )


test_that("calc_J works", {

  J <- calc_J(  d_m = "d2.2_m2rc",
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
    d_m = "d2.2_m2rc",
    typesample = "nbar",
    J = 5,
    MDES = 0.05, target.power = 0.80, tol = 0.01,
    Tbar = 0.50, alpha = 0.05, two.tailed = TRUE, numCovar.1 = 5,
    numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  )

  calcnbar <- pump_sample_raw( d_m = "d2.2_m2rc",
                               typesample = "nbar",
                               J = 10,
                               MDES = 0.05, target.power = 0.80,
                               Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                               numCovar.1 = 5, numCovar.2 = 1,
                               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcnbar
  expect_true( is.na( calcnbar$ss ) )


  calcJ <- pump_sample_raw( d_m = "d2.1_m2fc",
                            typesample = "J",
                            nbar = 258,
                            MDES = 0.05, target.power = 0.80,
                            Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                            numCovar.1 = 5, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcJ
  expect_true( !is.na( calcJ$ss ) )

  calcn <- pump_sample_raw( d_m = "d2.1_m2fc",
                            typesample = "nbar",
                            J = calcJ$ss,
                            MDES = 0.05, target.power = 0.80,
                            Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                            numCovar.1 = 5, numCovar.2 = 1,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcn
  expect_equal(258, calcn$ss, tol = 0.01)

  calcJ2 <- pump_sample_raw( d_m = "d2.1_m2fc",
                             typesample = "J",
                             nbar = calcn$ss,
                             MDES = 0.05, target.power = 0.80,
                             Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                             numCovar.1 = 5, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcJ2
  expect_equal(calcJ$ss, calcJ2$ss, tol = 0.01)

  calcn2 <- pump_sample_raw( d_m = "d2.1_m2fc",
                             typesample = "nbar",
                             J = calcJ2$ss,
                             MDES = 0.05, target.power = 0.80,
                             Tbar = 0.50, alpha = 0.05, two.tailed = TRUE,
                             numCovar.1 = 5, numCovar.2 = 1,
                             R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcn2
  expect_equal(calcn$ss, calcn2$ss, tol = 0.01)
})



test_that("BF for non individual power", {

  set.seed( 44941112 )
  p <- pump_power(  d_m = "d2.1_m2fc",
                    MTP = "BF",
                    J = 10,
                    nbar = 200,
                    M = 3,
                    MDES = rep(0.05, 3),
                    Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                    numCovar.1 = 5,
                    R2.1 = 0.1, ICC.2 = 0.05,
                    rho = 0.4, tnum = default.tnum )
  p

  ss <- pump_sample(    d_m = "d2.1_m2fc",
                        MTP = "BF",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "min1",
                        M = 3,
                        MDES = 0.05, target.power = p$min1[2],
                        tol = 0.01,
                        Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                        numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        rho = 0.4, tnum = default.tnum )


  expect_equal(ss$`Sample.size`, 10, tol = 1)
} )


test_that("plot_power_curve", {
    
  set.seed( 44941112 )
  ss1 <- pump_sample(   d_m = "d2.1_m2fc",
                        MTP = "BF",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        rho = 0.4 )
  expect_true(!is.null(plot_power_curve(ss1)))

  ss2 <- pump_sample(   d_m = "d2.1_m2fc",
                        MTP = "HO",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        rho = 0.4 )
  expect_true(!is.null(plot_power_curve(ss2)))
} )


test_that("plot_power_search", {
  ss1 <- pump_sample(   d_m = "d2.1_m2fc",
                        MTP = "BF",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        rho = 0.4 )
  expect_error(plot_power_search(ss1))

  ss2 <- pump_sample(   d_m = "d2.1_m2fc",
                        MTP = "HO",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        rho = 0.4 )
  expect_true(!is.null(plot_power_search(ss2)))
} )


test_that("pump_sample 2 level/2 level", {
  ss2 <- pump_sample(   d_m = "d2.1_m2fc",
                        MTP = "HO",
                        typesample = "J",
                        nbar = 200,
                        power.definition = "D1indiv",
                        M = 3,
                        MDES = 0.05, target.power = 0.80, tol = 0.01,
                        Tbar = 0.50, alpha = 0.05,
                        numCovar.1 = 5,
                        R2.1 = 0.1, ICC.2 = 0.05,
                        rho = 0.4 )
  ss2

  p2 <- pump_power( d_m = "d2.1_m2fc",
                    MTP = "HO",
                    J = ss2$`Sample.size`,
                    nbar = 200,
                    M = 3,
                    MDES = rep(0.05, 3),
                    Tbar = 0.50, alpha = 0.05,
                    numCovar.1 = 5,
                    R2.1 = 0.1,  ICC.2 = 0.05,
                    rho = 0.4 )

  p2
  expect_equal( p2[ 2, "indiv.mean" ], 0.80, tol = 0.02 )
} )


test_that("sample search when one end is missing", {

  set.seed( 20303 )
  pow_ref <- pump_power( d_m = "d2.2_m2rc",
                         MTP = "HO",
                         M = 4,
                         J = 10,
                         nbar = 10000,
                         MDES = rep( 0.40, 4 ),
                         Tbar = 0.50, alpha = 0.05,
                         numCovar.1 = 5, numCovar.2 = 1,
                         R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                         rho = 0.2 )

  pow_ref

  # this converges, but not to the correct value
  # because the power curve is too flat
  set.seed( 20303 )
  expect_warning( nbar1 <- pump_sample(
    d_m = "d2.2_m2rc",
    typesample = "nbar",
    power.definition = "min1",
    MTP = "HO",
    M = 4,
    J = 10,
    MDES = 0.40, target.power = pow_ref$min1[2], tol = 0.01,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
    rho = 0.2 ) )
  nbar1
  expect_true( !is.na( nbar1$`Sample.size` ) )


  # same problem happens with logit
  set.seed( 20303 )
  expect_warning( nbar2 <- pump_sample(
    d_m = "d2.2_m2rc",
    typesample = "nbar",
    power.definition = "min1",
    MTP = "HO",
    M = 4,
    J = 10,
    MDES = 0.40, target.power = pow_ref$min1[2], tol = 0.01,
    Tbar = 0.50, alpha = 0.05,
    numCovar.1 = 5, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
    rho = 0.2 ) )
  nbar2
  expect_true( !is.na( nbar2$`Sample.size` ) )

  # Now an infeasible calculation where the correlation makes min1 not able to
  # achieve power, even though independence would.
  set.seed( 443434344 )
  expect_warning(nbar3 <- pump_sample( d_m = "d2.2_m2rc",
                                          typesample = "nbar",
                                          power.definition = "min1",
                                          MTP = "HO",
                                          M = 4,
                                          J = 10,
                                          MDES = 0.39, target.power = 0.80, tol = 0.01,
                                          Tbar = 0.50, alpha = 0.05,
                                          numCovar.1 = 5, numCovar.2 = 1,
                                          R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                                          rho = 0.2, tnum = 1000 ) )
  nbar3
  expect_true( !is.na( nbar3$`Sample.size` ) )
  expect_true( nbar3$`Sample.size` > 100000 )
})


test_that("No adjustment", {

  nbar <- pump_sample(
    d_m = "d2.2_m2rc",
    MTP = 'BF',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = 0.8,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  )

  nbar <- expect_warning(pump_sample(
    d_m = "d2.2_m2rc",
    MTP = 'None',
    power.definition = 'D1indiv',
    typesample = 'nbar',
    target.power = 0.8,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  ))

  expect_error(nbar <- expect_warning(pump_sample(
    d_m = "d2.2_m2rc",
    MTP = 'None',
    power.definition = 'complete',
    typesample = 'nbar',
    target.power = 0.8,
    J = 60,
    M = 3,
    MDES = 0.125,
    Tbar = 0.5, alpha = 0.05, numCovar.1 = 1, numCovar.2 = 1,
    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, rho = 0.2
  )))

})

test_that( "sample errors out for MDES vector", {
  expect_error(pp <- pump_sample(
    d_m = "d2.2_m2rc",
    MTP = c("HO"),
    typesample = c("J"),
    MDES = rep(0.2, 5),
    M = 5,
    numZero = 1,
    nbar = 50,
    target.power = 0.8,
    power.definition = "indiv.mean",
    alpha = 0.5,
    Tbar = 0.8,
    numCovar.1 = 2,
    rho = 0.2))

  expect_error(pp <- pump_sample(
    d_m = "d2.2_m2rc",
    MTP = c("HO"),
    typesample = c("J"),
    M = 5,
    numZero = 1,
    nbar = 50,
    target.power = 0.8,
    power.definition = "indiv.mean",
    alpha = 0.5,
    Tbar = 0.8,
    numCovar.1 = 2,
    rho = 0.2))
})

test_that("M > 1 with MTP None", {
    
    ss <- expect_warning(pump_sample(
        d_m = "d2.1_m2fc",
        target.power = 0.8,
        power.definition = 'D1indiv',
        typesample = 'J',
        MTP = "None",
        MDES = 0.2,
        M = 3,
        nbar = 258,
        Tbar = 0.50, # prop Tx
        alpha = 0.05, # significance level
        numCovar.1 = 5,
        R2.1 = 0.1,
        ICC.2 = 0.05,
        rho = 0.4
    ))
    expect_true( nrow( ss ) == 1 )
    
    expect_error(expect_warning(pump_sample(
        d_m = "d2.1_m2fc",
        target.power = 0.8,
        power.definition = 'min2',
        typesample = 'J',
        MTP = "None",
        MDES = 0.2,
        M = 3,
        nbar = 258,
        Tbar = 0.50, # prop Tx
        alpha = 0.05, # significance level
        numCovar.1 = 5,
        R2.1 = 0.1,
        ICC.2 = 0.05,
        rho = 0.4
    )))
})

test_that("Sample with different correlations", {
    
    skip_on_cran()
    
    # zero correlation
    p <- pump_power(  d_m = "d2.1_m2fc",
                      MTP = "HO",
                      J = 10,
                      nbar = 200,
                      M = 20,
                      MDES = rep(0.05, 20),
                      Tbar = 0.50, alpha = 0.05,
                      numCovar.1 = 5,
                      R2.1 = 0.1, ICC.2 = 0.05,
                      rho = 0, tnum = default.tnum )
    p
    
    ss <- pump_sample(    d_m = "d2.1_m2fc",
                          MTP = "HO",
                          typesample = "J",
                          nbar = 200,
                          power.definition = "min1",
                          M = 20,
                          MDES = 0.05, target.power = p$min1[2],
                          tol = 0.01,
                          Tbar = 0.50, alpha = 0.05,
                          numCovar.1 = 5,
                          R2.1 = 0.1, ICC.2 = 0.05,
                          rho = 0 )
    
    
    expect_equal(ss$`Sample.size`, 10, tol = 1)
    
    
    # high correlation
    p <- pump_power(  d_m = "d2.1_m2fc",
                      MTP = "HO",
                      J = 10,
                      nbar = 200,
                      M = 20,
                      MDES = rep(0.05, 20),
                      Tbar = 0.50, alpha = 0.05,
                      numCovar.1 = 5,
                      R2.1 = 0.1, ICC.2 = 0.05,
                      rho = 0.95, tnum = default.tnum )
    p
    
    ss <- pump_sample(    d_m = "d2.1_m2fc",
                          MTP = "HO",
                          typesample = "J",
                          nbar = 200,
                          power.definition = "min1",
                          M = 20,
                          MDES = 0.05, target.power = p$min1[2],
                          tol = 0.01,
                          Tbar = 0.50, alpha = 0.05,
                          numCovar.1 = 5,
                          R2.1 = 0.1, ICC.2 = 0.05,
                          rho = 0.95 )
    
    
    expect_equal(ss$`Sample.size`, 10, tol = 1)
    
} )

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
        rho = 0.4, tnum = default.tnum )
    
    # sanity check: higher ICC means higher power
    expect_true(pow$D2indiv[1] > pow$D1indiv[1])
    expect_true(pow$D3indiv[1] > pow$D2indiv[1])
    
    ss1 <- pump_sample(
        d_m = "d2.1_m2fc",
        MTP = "HO",
        typesample = 'J',
        target.power = 0.8,
        power.definition = 'D1indiv',
        nbar = 200,
        M = 3,
        MDES = 0.05,
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 5,
        R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
        rho = 0.4
    )
    
    ss2 <- pump_sample(
        d_m = "d2.1_m2fc",
        MTP = "HO",
        typesample = 'J',
        target.power = 0.8,
        power.definition = 'D2indiv',
        nbar = 200,
        M = 3,
        MDES = 0.05,
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 5, 
        R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
        rho = 0.4
    )
    
    ss3 <- pump_sample(
        d_m = "d2.1_m2fc",
        MTP = "HO",
        typesample = 'J',
        target.power = 0.8,
        power.definition = 'D3indiv',
        nbar = 200,
        M = 3,
        MDES = 0.05,
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 5,
        R2.1 = 0.1, ICC.2 = c(0.1, 0.5, 0.8),
        rho = 0.4
    )
    
    # for same target power, we should need a bigger sample size for smaller ICC
    expect_true(ss1$`Sample.size` > ss2$`Sample.size`)
    expect_true(ss2$`Sample.size` > ss3$`Sample.size`)
    
})

test_that( "different MDES values work", {
    
    skip_on_cran()
    
    set.seed(034443)
    
    pow <- pump_power(
        d_m = "d2.2_m2rc",
        MTP = "BH",
        J = 40,
        nbar = 100,
        M = 5,
        MDES = c( 0.10, 0.05, 0.6, 0.05, 0.7 ),
        Tbar = 0.50, alpha = 0.05,
        numCovar.1 = 5, numCovar.2 = 1,
        R2.1 = 0.1, R2.2 = 0.7, ICC.2 = c(0.1, 0.5, 0.8, 0.01, 0.4),
        rho = 0.4, tnum = default.tnum ) 
    pow
    
    ss1 <- update( pow, type = "sample",
                   typesample = "J",
                   power.definition = "D4indiv",
                   target.power = pow$D4indiv[[2]],
                   max.steps = 30)
    ss1
    expect_equal( ss1$Sample.size, 40, tol = 0.05 )
    
    ss2 <- update( pow, type = "sample",
                   typesample = "J",
                   power.definition = "D5indiv",
                   target.power = 0.80 )
    ss2
    expect_true( ss2$Sample.size < 40 )
    
})
