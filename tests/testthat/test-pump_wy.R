set.seed(458738957)

source(here::here("R", "pump_wy.R"))

#### BASICS


########
# in this case the rawt values are too low
rawt <- c(2, 2, 2)
nullt <- c(3, 3, 3)

# what we want
exp.ind.B <- c(1, 1, 1)
# actual value
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
########


########
# in this case the rawt values are all high
rawt <- c(2, 2, 2)
nullt <- c(1, 1, 1)

# what we want
exp.ind.B <- c(0, 0, 0)
# actual value
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
########


########
# make sure the absolute value part is working
rawt <- c(-2, -2, -2)
nullt <- c(1, 1, 1)

# what we want
exp.ind.B <- c(0, 0, 0)
# actual value
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
########


########
# a mix of high and low values
rawt <- c(2, 0.5, 2)
nullt <- c(1, 1, 1)

# what we want
exp.ind.B <- c(0, 1, 0)
# actual value
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
########


#### test in-situ

# number of outcomes
M <- 3
# number of WY iterations
B <- 2
# degrees of freedom of null distribution
t.df <- 100

# correlation between outcomes
rho <- 0.2

# correlation matrix between outcomes
sigma <- matrix(rho, M, M)
diag(sigma) <- 1

# generate draws from null t
nullt <- mvtnorm::rmvt(B, sigma = sigma, df = t.df)

# > nullt
# [,1]       [,2]       [,3]
# [1,] -1.3717149 -0.9666078 -2.6004382
# [2,] -0.6419534 -0.4315923  0.9497158

# generate a rawt
rawt <- c(1.5, 2, 3)

exp.ind.B <- rbind(
  c(1, 1, 0),
  c(0, 0, 0)
)

ind.B <- t(apply(nullt, 1, comp.rawt.ss, rawt = rawt))

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

# now test p value
exp.adjp <- c(0.5, 0.5, 0)
adjp <- colMeans(ind.B)
test_that("p-values match", {
  expect_equal(adjp, exp.adjp)
})

## even more in situ

# generate a rawt
rawt.matrix <- rbind(
  c(1.5, 2, 3),
  c(0, -1.5, -3)
)

set.seed(458738957)
adjp <- adjp.wyss(rawt.matrix = rawt.matrix, B, sigma, t.df)

# null t
# t = 1
# [,1]       [,2]        [,3]
# [1,] -0.4223584 -0.4382459 -0.07559433
# [2,]  0.5050856 -2.6604902  1.17087901
# t = 2
# [,1]       [,2]      [,3]
# [1,] -0.5081732 -0.5777545 1.0057914
# [2,] -0.1043141 -2.8518170 0.1446717

exp.ind.B.1 <- rbind(
  c(0, 0, 0),
  c(1, 1, 0)
)
exp.pval.1 <- c(0.5, 0.5, 0)

exp.ind.B.2 <- rbind(
  c(1, 0, 0),
  c(1, 1, 0)
)
exp.pval.2 = c(1, 0.5, 0)

exp.adjp <- unname(rbind(exp.pval.1, exp.pval.2))
test_that("p-values match", {
  expect_equal(adjp, exp.adjp)
})

