# Return characteristics of a given context/d_m code (support function)

Returns number of levels and model at each level. See
pump_info()\$Context to get a list of supported d_ms.

## Usage

``` r
parse_d_m(d_m, d_only = FALSE)
```

## Arguments

- d_m:

  string; context to parse.

- d_only:

  TRUE/FALSE; TRUE means only look at design, ignore model if present.

## Value

list; list of features including number of levels, level of
randomization, etc.

## Examples

``` r
supported <- pump_info(comment = FALSE)$Context
parse_d_m( supported$d_m[4] )
#> $levels
#> [1] 2
#> 
#> $rand_level
#> [1] 1
#> 
#> $model2
#> [1] "fr"
#> 
#> $model2.p
#> [1] "f" "r"
#> 
#> $model3
#> NULL
#> 
#> $model3.p
#> NULL
#> 
#> $FE.2
#> [1] TRUE
#> 
#> $FE.3
#> [1] NA
#> 
#> $design
#> [1] "d2.1"
#> 
```
