# library( PUMP )
# library( testthat )

#--------------------------------------------------------------------------#
# single step
#--------------------------------------------------------------------------#


#----------------------------------------------------#
# test indicator function
#----------------------------------------------------#

#--------------------------#
rawt <- c(2, 2, 2)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(3, 3, 3)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)

exp.ind.B <- c(1, 1, 1)
ind.B <- comp_rawp_ss(nullp, rawp)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(2, 2, 2)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(1, 1, 1)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)

exp.ind.B <- c(0, 0, 0)
ind.B <- comp_rawp_ss(nullp, rawp)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#

#--------------------------#
rawt <- c(-2, -2, -2)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(1, 1, 1)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)

exp.ind.B <- c(0, 0, 0)
ind.B <- comp_rawp_ss(nullp, rawp)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(-2, -2, -2)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = FALSE)
nullt <- c(1, 1, 1)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = FALSE)

exp.ind.B <- c(1, 1, 1)
ind.B <- comp_rawp_ss(nullp, rawp)

test_that("Indicator functions match", {
    expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(2, 0.5, 2)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(1, 1, 1)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)

exp.ind.B <- c(0, 1, 0)
ind.B <- comp_rawp_ss(nullp, rawp)

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
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)

# > nullt
# [,1]       [,2]      [,3]
# [1,] -0.07034012 -0.2994425 0.1190887
# [2,]  0.83869418 -0.7146643 1.7726691

rawt <- c(0.2, 0.75, -2)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)

exp.ind.B <- rbind(
  c(1, 0, 0),
  c(1, 1, 0)
)
ind.B <- t(apply(nullp, 1, comp_rawp_ss, rawp = rawp))

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

rawt.mat <- rbind(
  c(1.5, 2, 3),
  c(0, -1.5, -3)
)
rawp.mat <- calc_pval(rawt.mat, t.df = 10, two.tailed = TRUE)

set.seed(458738957)
adjp <- adjp_wyss(rawp.mat = rawp.mat, B, sigma, t.df, two.tailed = TRUE)

# null p
# t = 1
# > nullp.mat
# [,1]        [,2]      [,3]
# [1,] 0.6736701 0.662153321 0.9398929
# [2,] 0.6146100 0.009090213 0.2444300
# t = 2
# > nullp.mat
# [,1]        [,2]      [,3]
# [1,] 0.6124511 0.564728575 0.3169435
# [2,] 0.9171290 0.005281889 0.8852613

exp.ind.B.1 <- rbind(
  c(0, 0, 0),
  c(1, 1, 1)
)
exp.pval.1 <- c(0.5, 0.5, 0.5)

exp.ind.B.2 <- rbind(
  c(1, 0, 0),
  c(1, 1, 1)
)
exp.pval.2 <- c(1, 0.5, 0.5)

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
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(3, 2, 4)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)
rawp.order <- order(rawp, decreasing = FALSE)

rawp.ordered <- rawp[rawp.order]
nullp.ordered <- nullp[rawp.order]

qstar <- rep(NA, M)
qstar[M] <- nullp.ordered[M]
for (h in (M-1):1)
{
    qstar[h] <- min(qstar[h + 1], nullp.ordered[h])
}

exp.qstar <- c(0.002518333, 0.013343655, 0.013343655)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(1, 1, 1)
ind.B <- comp_rawp_sd(nullp, rawp, rawp.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

#--------------------------#


#--------------------------#
rawt <- c(2.2, 2.1, 4.1)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(4, 3, 2)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)
rawp.order <- c(3, 1, 2)

rawp.ordered <- rawp[rawp.order]
nullp.ordered <- nullp[rawp.order]

qstar <- rep(NA, M)
qstar[M] <- nullp.ordered[M]
for (h in (M-1):1)
{
    qstar[h] <- min(qstar[h + 1], nullp.ordered[h])
}

exp.qstar <- c(0.002518333, 0.002518333, 0.013343655)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(0, 1, 1)
ind.B <- comp_rawp_sd(nullp, rawp, rawp.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

#--------------------------#


#--------------------------#
rawt <- c(-2.2, 2.1, 2.3)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(3, -2, 2.25)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)
rawp.order <- c(3, 1, 2)

rawp.ordered <- rawp[rawp.order]
nullp.ordered <- nullp[rawp.order]
qstar <- rep(NA, M)
qstar[M] <- nullp.ordered[M]
for (h in (M-1):1)
{
    qstar[h] <- min(qstar[h + 1], nullp.ordered[h])
}

exp.qstar <- c(0.01334366, 0.01334366, 0.07338803)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(1, 1, 0)
ind.B <- comp_rawp_sd(nullp, rawp, rawp.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})
#--------------------------#


#--------------------------#
rawt <- c(-3.3, 2.1, 2.3)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = FALSE)
nullt <- c(-3, -2, 2.25)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = FALSE)
rawp.order <- c(3, 2, 1)

rawp.ordered <- rawp[rawp.order]
nullp.ordered <- nullp[rawp.order]
qstar <- rep(NA, M)
qstar[M] <- nullp.ordered[M]
for (h in (M-1):1)
{
    qstar[h] <- min(qstar[h + 1], nullp.ordered[h])
}

exp.qstar <- c(0.02408983, 0.96330598, 0.99332817)
test_that("Qstar match", {
    expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(0, 0, 1)
ind.B <- comp_rawp_sd(nullp, rawp, rawp.order)

test_that("Indicator functions match", {
    expect_equal(ind.B, exp.ind.B)
})
#--------------------------#

#--------------------------#
rawt <- c(-5, 1.1, 4)
rawp <- calc_pval(rawt, t.df = 10, two.tailed = TRUE)
nullt <- c(3, -2, 1.8)
nullp <- calc_pval(nullt, t.df = 10, two.tailed = TRUE)
rawp.order <- c(1, 3, 2)

rawp.ordered <- rawp[rawp.order]
nullp.ordered <- nullp[rawp.order]
qstar <- rep(NA, M)
qstar[M] <- nullp.ordered[M]
for (h in (M-1):1)
{
    qstar[h] <- min(qstar[h + 1], nullp.ordered[h])
}

exp.qstar <- c(0.01334366, 0.07338803, 0.07338803)
test_that("Qstar match", {
  expect_equal(qstar, exp.qstar)
})

# test indicator B
exp.ind.B <- c(0, 0, 1)
ind.B <- comp_rawp_sd(nullp, rawp, rawp.order)

test_that("Indicator functions match", {
  expect_equal(ind.B, exp.ind.B)
})

#--------------------------#


#----------------------------------------------------#
# test monotonicity by itself
#----------------------------------------------------#

#--------------------------#
pi.p.m <- c(0.1, 0.2, 0.3)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp <- c(0.1, 0.2, 0.3)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#

#--------------------------#
pi.p.m <- c(0.1, 0.1, 0.2)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp <- c(0.1, 0.1, 0.2)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#


#--------------------------#
pi.p.m <- c(0.1, 0.2, 0.1)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp <- c(0.1, 0.2, 0.2)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#


#--------------------------#
pi.p.m <- c(0.1, 0.2, 0.1)

# enforce monotonicity (keep everything in same order as sorted RAW pvalues from original data)
adjp.minp <- rep(NA, length(pi.p.m))
adjp.minp[1] <- pi.p.m[1]
for (h in 2:length(pi.p.m)) {
  adjp.minp[h] <- max(pi.p.m[h], adjp.minp[h-1])
}

exp.adjp.minp <- c(0.1, 0.2, 0.2)

test_that("Montonicity is correct", {
  expect_equal(adjp.minp, exp.adjp.minp)
})
#--------------------------#


#----------------------------------------------------#
# test order matrix
#----------------------------------------------------#

rawt.mat <- rbind(
  c(0.8, 0.4, 3),
  c(0.3, 1, 0.1)
)
rawp.mat <- calc_pval(rawt.mat, t.df = 10, two.tailed = TRUE)

rawp.order.matrix <- t(apply(rawp.mat, 1, order, decreasing = FALSE))
exp.rawp.order.matrix <- rbind(
  c(3, 1, 2),
  c(2, 1, 3)
)
test_that("Order matrix matches", {
  expect_equal(rawp.order.matrix, exp.rawp.order.matrix)
})


#----------------------------------------------------#
# test final output is in correct order
#----------------------------------------------------#

#--------------------------#
rawp.order <- c(3, 1, 2)

rawp[rawp.order]
rawp[rawp.order][order(rawp.order)]

# make up a fake ind.B
ind.B <- rbind(
  c(0, 1, 0),
  c(0, 1, 1)
)

adjp <- get_adjp_minp(ind.B, rawp.order)
exp.adjp <- c(1, 1, 0)

test_that("P value output matches", {
  expect_equal(adjp, exp.adjp)
})

#----------------------------------------------------#
# test whole process
#----------------------------------------------------#

#--------------------------#
rawt.mat <- rbind(
  c(0.8, -0.4, 3),
  c(0.3, 1, 0.1)
)
rawp.mat <- calc_pval(rawt.mat, t.df = 10, two.tailed = TRUE)

B <- 100
set.seed(4335)
adjp.sd <- adjp_wysd(rawp.mat, B, sigma, t.df, two.tailed = TRUE, cl = NULL)
set.seed(4335)
adjp.ss <- adjp_wyss(rawp.mat, B, sigma, t.df, two.tailed = TRUE)
test_that("Smallest p-value should match", {
  expect_equal(adjp.sd[1,3], adjp.ss[1,3])
})
test_that("Smallest p-value should match", {
  expect_equal(adjp.sd[2,2], adjp.ss[2,2])
})
#--------------------------#

# calculate by hand
#--------------------------#
rawt.mat <- rbind(
  c(0.8, -0.4, 3),
  c(0.3, 1, 0.1)
)
rawp.mat <- calc_pval(rawt.mat, t.df = 10, two.tailed = TRUE)
B <- 2
# nullt
# [,1]       [,2]       [,3]
# [1,]  1.5400465 0.00675733 -0.1173385
# [2,] -0.3683377 0.21203830 -0.8790077
# [,1]          [,2]       [,3]
# [1,]  0.4936366 -0.0007528196 -0.4612104
# [2,] -1.0086968 -0.6327400480  0.8786953

set.seed(4334)
adjp <- adjp_wysd(rawp.mat, B, sigma, t.df, two.tailed = TRUE, cl = NULL)
exp.adjp <- rbind(
  c(0.5, 0.5, 0),
  c(1, 0.5, 1)
)

test_that("P value output matches", {
  expect_equal(adjp, exp.adjp)
})
#--------------------------#
