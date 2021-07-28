# library(testthat)
# source(here::here("R", "pump_power.R"))

#--------------------------------------------------------------------------#
# all outcomes are nonzero
#--------------------------------------------------------------------------#

pval.mat <- matrix(
  rbind(c(0,   0,   0  ),
        c(0,   0,   0.2),
        c(0,   0.2, 0.2),
        c(0.2, 0.2, 0.2)),
  nrow = 4, ncol = 3
)
ind.nonzero <- c(TRUE, TRUE, TRUE)
alpha <- 0.05

exp.power.results <- data.frame(
  D1indiv = 0.75,
  D2indiv = 0.5,
  D3indiv = 0.25,
  indiv.mean = 0.5,
  min1 = 0.75,
  min2 = 0.5,
  complete = 0.25
)
power.results.out <- get.power.results(pval.mat, ind.nonzero, alpha)

test_that("Power matches when all outcomes are nonzero", {
  expect_equal(exp.power.results, power.results.out)
})

#----------------------------------------------------#
# some outcomes are zero
#----------------------------------------------------#

pval.mat <- matrix(
  rbind(c(0,   0,     0,   0  ),
        c(0,   0.2,   0,   0.2),
        c(0,   0.2,   0,   0.2),
        c(0.2, 0.2,   0.2, 0.2)),
  nrow = 4, ncol = 4
)
ind.nonzero <- c(TRUE, TRUE, FALSE, FALSE)
alpha <- 0.05

exp.power.results <- data.frame(
  D1indiv = 0.75,
  D2indiv = 0.25,
  indiv.mean = 0.5,
  min1 = 0.75,
  complete = 0.25
)
power.results.out <- get.power.results(pval.mat, ind.nonzero, alpha)

test_that("Power matches with a mix of zero and nonzero outcomes", {
  expect_equal(exp.power.results, power.results.out)
})
