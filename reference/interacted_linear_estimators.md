# Interacted linear regression models

Code taken from:
https://github.com/lmiratrix/blkvar/blob_master/R/linear_model_method.R

## Usage

``` r
interacted_linear_estimators(
  Yobs,
  Z,
  B,
  siteID = NULL,
  data = NULL,
  control_formula = NULL,
  use.lmer = FALSE
)
```

## Value

Dataframe of the different versions of this estimator (person and site
weighted)

## Details

These linear models have block by treatment interaction terms. The final
ATE estimates are then weighted average of the block (site) specific ATE
estimates.

If siteID passed, it will weight the RA blocks within site and then
average these site estimates.

SEs come from the overall variance-covariance matrix.
