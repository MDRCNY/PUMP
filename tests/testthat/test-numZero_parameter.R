# library( PUMP )
# library( testthat )


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# ----- three level models ------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# --------    d3.1_m3rr2rr    --------
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #

skip_on_cran()

test_that("numZero reduces power as expected and can be inverted", {

  set.seed(101003)

  pp <- pump_power( d_m = "d2.2_m2rc",
                   MTP = "HO",
                   M = 4,
                   J = 30,
                   nbar = 100,
                   MDES = rep( 0.15, 4),
                   Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                   numCovar.1 = 5, numCovar.2 = 1,
                   R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                   rho = 0.2, tnum=10000 )
  pp
  
  pp2 <- pump_power( d_m = "d2.2_m2rc",
                    MTP = "HO",
                    M = 14,
                    numZero = 10,
                    J = 30,
                    nbar = 100,
                    MDES = rep( 0.15, 4),
                    Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                    numCovar.1 = 5, numCovar.2 = 1,
                    R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                    rho = 0.2, tnum=10000, verbose = FALSE )
  pp2
  
  expect_equal( as.data.frame(pp)[1,2:5], 
                as.data.frame(pp2)[1,2:5], tolerance = 0.05 )
  
  expect_true( all( as.data.frame(pp)[2,2:10] >= as.data.frame(pp2)[2,2:10] ) ) 
  

  pp2flipHand <- pump_mdes( d_m = "d2.2_m2rc",
                           MTP = "HO",
                           target.power = pp2$min2[[2]], power.definition = "min2",
                           M = 14,
                           numZero = 10,
                           J = 30,
                           nbar = 100,
                           Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                           numCovar.1 = 5, numCovar.2 = 1,
                           R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                           rho = 0.2 )
  pp2flipHand
  
  pp2flip <- update( pp2, type="mdes", target.power = pp2$min2[[2]], power.definition = "min2" )
  pp2flip
  
  expect_equal( pp2flipHand$Adjusted.MDES,  pp2flip$Adjusted.MDES, tolerance = 0.01 )
  expect_equal( pp2flip$Adjusted.MDES, 0.15, tolerance = 0.01 )
  

  
  #### Sample size flip #####
  
  
  pp2

  ppSS <- pump_sample( d_m = "d2.2_m2rc",
                      MTP = "HO",
                      typesample = "J",
                      power.definition = "min2",
                      target.power = pp2$min2[[2]],
                      M = 14,
                      numZero = 10,
                      MDES = 0.15,
                      nbar = 100,
                      Tbar = 0.50, alpha = 0.05, two.tailed = FALSE,
                      numCovar.1 = 5, numCovar.2 = 1,
                      R2.1 = 0.1, R2.2 = 0.7, ICC.2 = 0.05,
                      rho = 0.2, verbose = FALSE )
  ppSS
  
  expect_equal( ppSS$Sample.size, 30, tolerance = 1 )

} )
