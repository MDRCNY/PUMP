# Fit a bounded logistic curve

Curve is of form f(y) = pmin + (pmax-pmin) \* logistic( beta0 + beta1\*x
)

## Usage

``` r
fit_bounded_logistic(x, y, wt)
```

## Arguments

- x:

  The vector of covariate values of the logistics

- y:

  The proportion of 1s for a given value of x. Same length as x.

- wt:

  The weight to place on a given x-y pair. Same length as x, or scalar.

## Value

Vector of four estimated parameters for the logistic curve: beta0,
beta1, pmin, pmax

## Details

(logistic as defined by plogis) Note that a logistic curve is not a
perfect fit for the functional form of the power curve, but is a useful
approximation for the search procedure.
