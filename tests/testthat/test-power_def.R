# library( PUMP )
# library( testthat )

#--------------------------------------------------------------------------#
# all outcomes are nonzero
#--------------------------------------------------------------------------#

unadj.pval.mat <- adj.pval.mat <- matrix(
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
power.results.out <- get_power_results(adj.pval.mat, unadj.pval.mat,
                                       ind.nonzero, alpha)

test_that("Power matches when all outcomes are nonzero", {
  expect_equal(power.results.out, exp.power.results)
})

#----------------------------------------------------#
# some outcomes are zero
#----------------------------------------------------#

unadj.pval.mat <- adj.pval.mat  <- matrix(
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
  min2 = 0.75,
  complete = NA
)
power.results.out <- get_power_results(adj.pval.mat, unadj.pval.mat,
                                       ind.nonzero, alpha)

test_that("Power matches with a mix of zero and nonzero outcomes", {
  expect_true(is.na(power.results.out$complete))
  expect_equal(power.results.out[,1:5], exp.power.results[,1:5])
})


#----------------------------------------------------#
# some outcomes are zero part II
#----------------------------------------------------#

unadj.pval.mat <- adj.pval.mat  <- matrix(
  rbind(c(0,     0,     0,     0  ),
        c(0.1,   0.2,   0,     0.2),
        c(0,     0.2,   0,     0.2),
        c(0.2,   0.2,   0.2,   0.2)),
  nrow = 4, ncol = 4
)
ind.nonzero <- c(TRUE, TRUE, TRUE, FALSE)
alpha <- 0.05

exp.power.results <- data.frame(
  D1indiv = 0.5,
  D2indiv = 0.25,
  D3indiv = 0.75,
  indiv.mean = 0.5,
  min1 = 0.75,
  min2 = 0.5,
  min3 = 0.25,
  complete = NA
)
power.results.out <- get_power_results(adj.pval.mat, unadj.pval.mat,
                                       ind.nonzero, alpha)

test_that("Power matches with a mix of zero and nonzero outcomes", {
  expect_true(is.na(power.results.out$complete))
  expect_equal(power.results.out[,1:6], exp.power.results[,1:6])
})


#----------------------------------------------------#
# complete power with unadjusted
#----------------------------------------------------#

unadj.pval.mat <- matrix(
  rbind(c(0,     0,     0,     0  ),
        c(0.1,   0.2,   0,     0.2),
        c(0,     0,     0,     0.2),
        c(0.2,   0.2,   0.2,   0.2)),
  nrow = 4, ncol = 4
)
adj.pval.mat <- matrix(
  rbind(c(0,     0,     0,     0  ),
        c(0.1,   0.2,   0.1,   0.2),
        c(0,     0.2,   0,     0.2),
        c(0.2,   0.2,   0.2,   0.2)),
  nrow = 4, ncol = 4
)

ind.nonzero <- c(TRUE, TRUE, TRUE, TRUE)
alpha <- 0.05

exp.power.results <- data.frame(
  D1indiv = 0.5,
  D2indiv = 0.25,
  D3indiv = 0.5,
  D4indiv = 0.25,
  indiv.mean = 0.375,
  min1 = 0.5,
  min2 = 0.5,
  min3 = 0.25,
  complete = 0.25
)
power.results.out <- get_power_results(adj.pval.mat, unadj.pval.mat,
                                       ind.nonzero, alpha)

test_that("Power matches with a mix of zero and nonzero outcomes", {
  expect_equal(power.results.out, exp.power.results)
})


#----------------------------------------------------#
# perfectly correlated outcomes
#----------------------------------------------------#

unadj.pval.mat <- adj.pval.mat <- matrix(
    rbind(c(0,    0,     0,    0  ),
          c(0.2,  0.2,   0.2,  0.2),
          c(0.2,  0.2,   0.2,  0.2),
          c(0.2,  0.2,   0.2,  0.2)),
    nrow = 4, ncol = 4
)
ind.nonzero <- c(TRUE, TRUE, TRUE, TRUE)
alpha <- 0.05

exp.power.results <- data.frame(
    D1indiv = 0.25,
    D2indiv = 0.25,
    D3indiv = 0.25,
    D4indiv = 0.25,
    indiv.mean = 0.25,
    min1 = 0.25,
    min2 = 0.25,
    min3 = 0.25,
    complete = 0.25
)
power.results.out <- get_power_results(adj.pval.mat, unadj.pval.mat,
                                       ind.nonzero, alpha)

test_that("Power matches with perfectly correlated outcomes", {
    expect_equal(power.results.out, exp.power.results)
})
