R2.1 = 0.61
pmdes <- pump_mdes(design = "d1.1_m1c", MTP = "Holm",
                   target.power = 0.80, power.definition = "min1", tol = 0.02,
                   R2.1 = R2.1, numCovar.1 = 1, J = 1,
                   max.tnum = 1000,
                   M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5)
pmdes



a = pump_mdes(design = "d1.1_m1c", MTP = "Bonferroni",
          target.power = 0.80,
          power.definition = "min1", tol = 0.02,
          R2.1 = R2.1, numCovar.1 = 1, J = 1,
          max.tnum = 1000,
          M = 3, nbar = 12, Tbar = 1/3, alpha = 0.10, rho = 0.5)
a

update( a, type = "sample", typesample = "nbar", MDES = 0.20 )

xx = update( a, type = "power", MDES = 0.20 )
xx
str(xx)

pp <- pump_mdes_grid(    design = "d3.2_m3ff2rc",
                         MTP = c( "Bonferroni", "Holm" ), 
                         target.power = c( 0.50, 0.80 ),
                         power.definition = c( "min1", "D1indiv" ),
                         tol = 0.05,
                         M = 5,
                         J = 5,
                         K = 7, # number RA blocks
                         nbar = 58,
                         Tbar = 0.50, # prop Tx
                         alpha = 0.15, # significance level
                         numCovar.1 = 1, numCovar.2 = 1,
                         R2.1 = 0.1, R2.2 = 0.7,
                         ICC.2 = 0.05, ICC.3 = 0.9,
                         rho = 0.4, # how correlated outcomes are
                         verbose = FALSE, max.tnum = 500,
)
pp



pp <- pump_sample_grid(    design = "d3.2_m3ff2rc",
                         MTP = c( "Bonferroni", "Holm" ), 
                         target.power = c( 0.50, 0.80 ),
                         typesample = "nbar",
                         MDES = 0.20,
                         power.definition = c( "min1", "D1indiv" ),
                         tol = 0.05,
                         M = 5,
                         J = 5,
                         K = 7, # number RA blocks
                         Tbar = 0.50, # prop Tx
                         alpha = 0.15, # significance level
                         numCovar.1 = 1, numCovar.2 = 1,
                         R2.1 = 0.1, R2.2 = 0.7,
                         ICC.2 = 0.05, ICC.3 = 0.9,
                         rho = 0.4, # how correlated outcomes are
                         verbose = FALSE, max.tnum = 500,
)
pp
