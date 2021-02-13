# library(testthat)
# source(here::here("R", "pump_wy.R"))

#--------------------------------------------------------------------------#
# single step
#--------------------------------------------------------------------------#


#----------------------------------------------------#
# test indicator function
#----------------------------------------------------#

#--------------------------#
rawt <- c(2, 2, 2)
nullt <- c(3, 3, 3)

exp.ind.B <- c(1, 1, 1)
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(2, 2, 2)
nullt <- c(1, 1, 1)

exp.ind.B <- c(0, 0, 0)
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(-2, -2, -2)
nullt <- c(1, 1, 1)

exp.ind.B <- c(0, 0, 0)
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(2, 0.5, 2)
nullt <- c(1, 1, 1)

exp.ind.B <- c(0, 1, 0)
ind.B <- comp.rawt.ss(nullt, rawt)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#----------------------------------------------------#
# test indicator function
#----------------------------------------------------#

#--------------------------#

M <- 3
B <- 2
t.df <- 100
rho <- 0.2
sigma <- matrix(rho, M, M)
diag(sigma) <- 1

# generate draws from null t
set.seed(433)
nullt <- mvtnorm::rmvt(B, sigma = sigma, df = t.df)

# > nullt
# [,1]       [,2]      [,3]
# [1,] -0.07034012 -0.2994425 0.1190887
# [2,]  0.83869418 -0.7146643 1.7726691

rawt <- c(0.2, 0.75, -2)

exp.ind.B <- rbind(
  c(1, 0, 0),
  c(1, 1, 0)
)
ind.B <- t(apply(nullt, 1, comp.rawt.ss, rawt = rawt))

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

exp.adjp <- c(1, 0.5, 0)
adjp <- colMeans(ind.B)
test_that("p-values match", {
  expect_equal(adjp, exp.adjp)
})

#----------------------------------------------------#
# test p values across multiple iterations
#----------------------------------------------------#

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

#--------------------------------------------------------------------------#
# step down
#--------------------------------------------------------------------------#


#----------------------------------------------------#
# test qstar and indB
#----------------------------------------------------#

#--------------------------#
rawt <- c(2.1, 2.2, 2.3)
nullt <- c(3, 2, 4)
rawt.order <- order(abs(rawt))

rawt.ordered <- rawt[rawt.order]
nullt.ordered <- abs(nullt[rawt.order])
qstar <- rep(NA, M)
qstar[1] <- nullt.ordered[1]
for (h in 2:M)
{
  qstar[h] <- max(abs(qstar[h - 1]), nullt.ordered[h])
}

exp.qstar <- c(3, 3, 4)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(1, 1, 1)
ind.B <- comp.rawt.sd(nullt, rawt, rawt.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

#--------------------------#


#--------------------------#
rawt <- c(2.2, 2.1, 2.3)
nullt <- c(3, 2, 4)
rawt.order <- c(2, 1, 3)

rawt.ordered <- rawt[rawt.order]
nullt.ordered <- abs(nullt[rawt.order])
qstar <- rep(NA, M)
qstar[1] <- nullt.ordered[1]
for (h in 2:M)
{
  qstar[h] <- max(abs(qstar[h - 1]), nullt.ordered[h])
}

exp.qstar <- c(2, 3, 4)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(0, 1, 1)
ind.B <- comp.rawt.sd(nullt, rawt, rawt.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

#--------------------------#


#--------------------------#
rawt <- c(-2.2, 2.1, 2.3)
nullt <- c(3, -2, 4)
rawt.order <- c(2, 1, 3)

rawt.ordered <- rawt[rawt.order]
nullt.ordered <- abs(nullt[rawt.order])
qstar <- rep(NA, M)
qstar[1] <- nullt.ordered[1]
for (h in 2:M)
{
  qstar[h] <- max(abs(qstar[h - 1]), nullt.ordered[h])
}
exp.qstar <- c(2, 3, 4)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(0, 1, 1)
ind.B <- comp.rawt.sd(nullt, rawt, rawt.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(-5, 1.1, 4)
nullt <- c(3, -2, 1.8)
rawt.order <- c(2, 3, 1)

rawt.ordered <- rawt[rawt.order]
nullt.ordered <- abs(nullt[rawt.order])
qstar <- rep(NA, M)
qstar[1] <- nullt.ordered[1]
for (h in 2:M)
{
  qstar[h] <- max(abs(qstar[h - 1]), nullt.ordered[h])
}
exp.qstar <- c(2, 2, 3)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(1, 0, 0)
ind.B <- comp.rawt.sd(nullt, rawt, rawt.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

#--------------------------#


#----------------------------------------------------#
# test monotonicity by itself
#----------------------------------------------------#

#--------------------------#
pi.p.m = c(0.1, 0.2, 0.3)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp = c(0.1, 0.2, 0.3)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#

#--------------------------#
pi.p.m = c(0.1, 0.1, 0.2)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp = c(0.1, 0.1, 0.2)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#


#--------------------------#
pi.p.m = c(0.1, 0.2, 0.1)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp = c(0.1, 0.2, 0.2)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#


#--------------------------#
pi.p.m = c(0.1, 0.2, 0.1)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp = c(0.1, 0.2, 0.2)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#


#----------------------------------------------------#
# test order matrix
#----------------------------------------------------#

rawt.matrix = rbind(
  c(0.8, 0.4, 3),
  c(0.3, 1, 0.1)
)

rawt.order.matrix <- t(apply(abs(rawt.matrix), 1, order, decreasing = FALSE))
exp.rawt.order.matrix <- rbind(
  c(2, 1, 3),
  c(3, 1, 2)
)
test_that("Order matrix matches", {
  expect_equal(rawt.order.matrix, exp.rawt.order.matrix)
})


#----------------------------------------------------#
# test final output is in correct order
#----------------------------------------------------#

#--------------------------#
rawt.order <- c(3, 1, 2)

rawt[rawt.order]
rawt[rawt.order][order(rawt.order)]

# make up a fake ind.B
ind.B <- rbind(
  c(0, 1, 0),
  c(0, 1, 1)
)

adjp <- get.adjp.minp(ind.B, rawt.order)
exp.adjp <- c(1, 1, 0)

test_that("P value output matches", {
  expect_equal(adjp, exp.adjp)
})


