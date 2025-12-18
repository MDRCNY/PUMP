# generate a parameterized covariance matrix from the provided 3 blocks

generate a parameterized covariance matrix from the provided 3 blocks

## Usage

``` r
gen_RE_cov_matrix(Sigma.w, Sigma.z, Sigma.wz)
```

## Arguments

- Sigma.w:

  level 3 covariance matrix

- Sigma.z:

  level 2 covariance matrix

- Sigma.wz:

  covariance between level 2 and level 3

## Value

2M x 2M matrix for generating correlated pairs of random effects
