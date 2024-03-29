% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correlation.R
\name{check_cor}
\alias{check_cor}
\title{Check correlation of test statistics (simulation function)}
\usage{
check_cor(
  pump.object = NULL,
  rho.V = NULL,
  rho.w0 = NULL,
  rho.w1 = NULL,
  rho.X = NULL,
  rho.u0 = NULL,
  rho.u1 = NULL,
  rho.C = NULL,
  rho.r = NULL,
  d_m = NULL,
  model.params.list = NULL,
  Tbar = 0.5,
  n.sims = 100
)
}
\arguments{
\item{pump.object}{A pumpresult object.}

\item{rho.V}{matrix; correlation matrix of level 3 covariates.}

\item{rho.w0}{matrix; correlation matrix of level 3 random effects.}

\item{rho.w1}{matrix; correlation matrix of level 3 random impacts.}

\item{rho.X}{matrix; correlation matrix of level 2 covariates.}

\item{rho.u0}{matrix; correlation matrix of level 2 random effects.}

\item{rho.u1}{matrix; correlation matrix of level 2 random impacts.}

\item{rho.C}{matrix; correlation matrix of level 1 covariates.}

\item{rho.r}{matrix; correlation matrix of level 1 residuals.}

\item{d_m}{string; a single context, which is a design and model
code. See pump_info() for list of choices.}

\item{model.params.list}{list; model parameters such as ICC, R2, etc. See
simulation vignette for details.}

\item{Tbar}{scalar; the proportion of samples that are assigned to
the treatment.}

\item{n.sims}{numeric; Number of simulated datasets to generate.
More datasets will achieve a more accurate result
but also increase computation time.}
}
\value{
matrix; M x M correlation matrix between test statistics.
}
\description{
Estimates the pairwise correlations
between test statistics for all outcomes.

Takes in two options:
- a pumpresult object
OR
- a list of necessary data-generating parameters
- the context (d_m)
- Tbar

Note that this function can take several minutes to run.
}
\examples{
pp <- pump_power( d_m = "d3.2_m3ff2rc",
                  MTP = "BF",
                  MDES = rep( 0.10, 2 ),
                  M = 2,
                  J = 4, # number of schools/block
                  K = 10, # number RA blocks
                  nbar = 50,
                  Tbar = 0.50, # prop Tx
                  alpha = 0.05, # significance level
                  numCovar.1 = 5, numCovar.2 = 3,
                  R2.1 = 0.1, R2.2 = 0.7,
                  ICC.2 = 0.05, ICC.3 = 0.4,
                  rho = 0.4, # how correlated test statistics are
                  tnum = 200
)
cor.tstat <- check_cor(
    pump.object = pp, n.sims = 4
)
est.cor <- mean(cor.tstat[lower.tri(cor.tstat)])
}
