# Helper function for Westfall Young Single Step

Used for the Westfall-Young single-step multiple testing procedure
(MTP). It compares whether any of the null values across outcomes
exceeds each raw value for each outcome

## Usage

``` r
comp_rawp_ss(nullp, rawp)
```

## Arguments

- nullp:

  a vector of p values under H0

- rawp:

  a vector of raw p values under H1

## Value

returns a vector of 1s and 0s with length of M outcomes
