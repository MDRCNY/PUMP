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

# calc rawp
rawp <- apply(X, 1, function(x){ return(t.test(x[T.x == 1], x[T.x == 0])$p.value) })

# multtest
set.seed(0217)
# suppress output
sink("/dev/null")
mult.out <- multtest::mt.maxT(X, classlabel = T.x, B = B)
mult.adjp <- mult.out$adjp[order(mult.out$index)]
sink()

# my version
rawp.order <- order(rawp, decreasing = FALSE)
set.seed(0217)
# get nullt
nullp <- NULL
# get nullt
for(b in 1:B)
{
  T.x.b <- sample(c(rep(1, 50), rep(0, 50)))
  nullp.b <- apply(X, 1, function(x){ return(t.test(x[T.x.b == 1], x[T.x.b == 0])$p.value) })
  nullp <- rbind(nullp, nullp.b)
}

ind.B <- t(apply(nullp, 1, comp_rawp_sd, rawp = rawp, rawp.order = rawp.order))
exp.adjp <- get_adjp_minp(ind.B, rawp.order)

test_that("P value output matches within 1%", {
  expect_equal(mult.adjp, exp.adjp, tol = 0.01)
})
