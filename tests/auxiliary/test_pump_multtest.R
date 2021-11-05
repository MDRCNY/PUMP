#----------------------------------------------------#
# compare to multtest
#----------------------------------------------------#

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("multtest")
library(multtest)

set.seed(0217)
B <- 10000

# generate fake input data
M <- 3
N <- 100
# treatment vector
T.x <- sample(c(rep(1, 50), rep(0, 50)))
# X has a treatment effect
X <- matrix(NA, nrow = M, ncol = N)
X[1,] <- 2 * T.x + rnorm(N)
X[2,] <- 0 * T.x + rnorm(N)
X[3,] <- -4 * T.x + rnorm(N)

# calc rawt
rawt <- apply(X, 1, function(x){ return(t.test(x[T.x == 1], x[T.x == 0])$statistic) })


# multtest
set.seed(0217)
# suppress output
sink("/dev/null")
mult.out <- multtest::mt.maxT(X, classlabel = T.x, B = B)
mult.adjp <- mult.out$adjp[order(mult.out$index)]
sink()

# my version
rawt.order <- order(abs(rawt), decreasing = TRUE)
set.seed(0217)
# get nullt
nullt <- NULL
# get nullt
for(b in 1:B)
{
  T.x.b <- sample(c(rep(1, 50), rep(0, 50)))
  nullt.b <- apply(X, 1, function(x){ return(t.test(x[T.x.b == 1], x[T.x.b == 0])$statistic) })
  nullt <- rbind(nullt, nullt.b)
}

ind.B <- t(apply(nullt, 1, comp.rawt.sd, rawt = rawt, rawt.order = rawt.order))
exp.adjp <- get.adjp.minp(ind.B, rawt.order)

test_that("P value output matches within 5%", {
  expect_equal(mult.adjp, exp.adjp, tol = 0.05)
})
