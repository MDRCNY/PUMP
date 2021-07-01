
library( pum )
library( testthat )


test_that("calc.nbar works", {

  nbar <- pum:::calc.nbar( design="simple_c2_2r",
                     MT = 2.8,
                     MDES = 0.20,
                     J = 5,
                     Tbar = 0.25,
                     R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  nbar
  expect_true( is.na( nbar ) )

  nbar <- pum:::calc.nbar( design="simple_c2_2r",
                     MT = 2.8,
                     MDES = 0.20,
                     J = 305,
                     Tbar = 0.25,
                     R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  nbar
  expect_true( !is.na( nbar ) )

} )




test_that("calc.J works", {

  J <- pum:::calc.J( design="simple_c2_2r",
          MT = 2.8,
                   MDES = 0.20,
                   nbar = 200,
                   Tbar = 0.25,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )
  J
  expect_true( J < 100 )
} )



test_that("pump_sample_raw works", {

  expect_error( calcnbar <- pump_sample_raw( design="simple_c2_2r",
                            typesample = "nbar",
                            J = 5,
                            MDES = 0.05, target.power = 0.80, tol = 0.01,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  )

  calcnbar <- pump_sample_raw( design="simple_c2_2r",
                               typesample = "nbar",
                               J = 10,
                               MDES = 0.05, target.power = 0.80, tol = 0.01,
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcnbar
  expect_true( is.na( calcnbar ) )


  calcJ <- pump_sample_raw( design="blocked_i1_2c",
                            typesample = "J",
                            nbar = 258,
                            MDES = 0.05, target.power = 0.80, tol = 0.01,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcJ
  expect_true( !is.na( calcJ ) )

  calcn <- pump_sample_raw( design="blocked_i1_2c",
                               typesample = "nbar",
                               J = calcJ,
                               MDES = 0.05, target.power = 0.80, tol = 0.01,
                               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05 )

  calcn
  expect_true( abs( calcn - 245 ) < 10 )

  calcJ2 <- pump_sample_raw( design="blocked_i1_2c",
                            typesample = "J",
                            nbar = calcn,
                            MDES = 0.05, target.power = 0.80, tol = 0.01,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcJ2
  # expect_true( abs( calcn - 258 ) < 10 )

  calcn2 <- pump_sample_raw( design="blocked_i1_2c",
                            typesample = "nbar",
                            J = calcJ2,
                            MDES = 0.05, target.power = 0.80, tol = 0.01,
                            Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
                            R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4 )

  calcn2
  #expect_true( abs( calcn - calcn2 ) < 10 )

  # THERE IS A PROBLEM HERE!  Shouldnt these two things align???
})


test_that("pump_sample 2 level/2 level", {
  pump_sample( design="blocked_i1_2c",
               MTP = "Holm",
               typesample = "J",
               nbar = 200,
               power.definition = "D1indiv",
               M = 3,
               MDES = 0.05, target.power = 0.80, tol = 0.01,
               Tbar = 0.50, alpha = 0.05, numCovar.1 = 5, numCovar.2 = 0,
               R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05, ICC.3 = 0.4,
               rho = 0.4 )

  expect_equal( dim( pp ), c(2,7) )

  expect_true( all( pp[,"min1"] >= pp[,"D1indiv"] ) )
})




