# pumpresult object for results of power calculations

The pumpresult object is an S3 class that holds the results from
\`pump_power()\`, \`pump_sample()\`, and \`pump_mdes()\`.

It has several methods that pull different information from this object,
and some printing methods for getting nicely formatted results.

Pump result objects are also data.frames, so they can be easily
manipulated and combined. The return values from the \`grid\` functions
will just return data frames in general.

Returns whether call was power, mdes, or sample.

Calls the print_context method with results and control both set to
TRUE.

## Usage

``` r
params(x, ...)

d_m(x, ...)

design(x, ...)

search_path(x, ...)

pump_type(x)

is.pumpresult(x)

# S3 method for class 'pumpresult'
x[...]

# S3 method for class 'pumpresult'
x[[...]]

# S3 method for class 'pumpresult'
dim(x, ...)

# S3 method for class 'pumpresult'
summary(object, ...)

# S3 method for class 'pumpresult'
print(x, n = 10, header = TRUE, search = FALSE, include_SE = TRUE, ...)

# S3 method for class 'pumpresult'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  a pumpresult object (except for is.pumpresult, where it is a generic
  object to check).

- ...:

  additional arguments to be passed to the as.data.frame.list methods.

- object:

  Object to summarize.

- n:

  Number of lines of search path to print, max.

- header:

  FALSE means skip some header info on the result, just print the
  data.frame of actual results.

- search:

  FALSE means don't print the search path for a result for mdes or
  sample.

- include_SE:

  TRUE means include standard errors given design (if any) in the
  printout. Default to TRUE.

- row.names:

  NULL or a character vector giving the row names for the data frame.

- optional:

  logical. If TRUE, setting row names and converting column names is
  optional.

## Value

params: List of design parameters used.

d_m: Context (d_m) used (as string).

design (the randomization and levels) as string.

search_path: Dataframe describing search path, if it was saved in the
pumpresult object.

pump_type: power, mdes, or sample, as a string.

is.pumpresult: TRUE if object is a pumpresult object.

\`\[\`: pull out rows and columns of the dataframe.

\`\[\[\`: pull out single element of dataframe.

dim: Dimension of pumpresult (as matrix)

summary: No return value; prints results.

print: No return value; prints results.

as.data.frame: pumpresult object as a clean dataframe (no more
attributes from pumpresult).

## See also

update

update_grid

print_context

print_context

## Examples

``` r
pp <- pump_power(d_m = "d3.2_m3ff2rc",
  MTP = 'HO', nbar = 50, J = 30, K = 10,
  M = 5, MDES = 0.125, Tbar = 0.5, alpha = 0.05,
  numCovar.1 = 1, numCovar.2 = 1,
  R2.1 = 0.1, R2.2 = 0.1, ICC.2 = 0.2, ICC.3 = 0.2,
  omega.2 = 0, omega.3 = 0.1, rho = 0.5, tnum = 1000)
  
print(pp)
#> power result: d3.2_m3ff2rc d_m with 5 outcomes
#>   MTP    D1indiv    D2indiv    D3indiv    D4indiv    D5indiv indiv.mean  min1
#>  None      0.670      0.692      0.683      0.681      0.712      0.688      
#>    SE ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )                  
#>    HO      0.527      0.544      0.541      0.541      0.557      0.542 0.784
#>   min2  min3  min4 complete df1
#>                             279
#>                                
#>  0.651 0.519 0.423    0.358    
#>  0.005 <= MCSE <= 0.008
params(pp)
#> $MTP
#> [1] "HO"
#> 
#> $MDES
#> [1] 0.125 0.125 0.125 0.125 0.125
#> 
#> $numZero
#> NULL
#> 
#> $M
#> [1] 5
#> 
#> $J
#> [1] 30
#> 
#> $K
#> [1] 10
#> 
#> $nbar
#> [1] 50
#> 
#> $Tbar
#> [1] 0.5
#> 
#> $alpha
#> [1] 0.05
#> 
#> $two.tailed
#> [1] TRUE
#> 
#> $numCovar.1
#> [1] 1
#> 
#> $numCovar.2
#> [1] 1
#> 
#> $numCovar.3
#> [1] 0
#> 
#> $R2.1
#> [1] 0.1 0.1 0.1 0.1 0.1
#> 
#> $R2.2
#> [1] 0.1 0.1 0.1 0.1 0.1
#> 
#> $R2.3
#> [1] 0 0 0 0 0
#> 
#> $ICC.2
#> [1] 0.2 0.2 0.2 0.2 0.2
#> 
#> $ICC.3
#> [1] 0.2 0.2 0.2 0.2 0.2
#> 
#> $omega.2
#> [1] 0 0 0 0 0
#> 
#> $omega.3
#> [1] 0.1 0.1 0.1 0.1 0.1
#> 
#> $rho
#> [1] 0.5
#> 
#> $rho.matrix
#> NULL
#> 
#> $B
#> [1] 1000
#> 
#> $tnum
#> [1] 1000
#> 
#> $d_m
#> [1] "d3.2_m3ff2rc"
#> 
print_context(pp)
#> power result: d3.2_m3ff2rc d_m with 5 outcomes
#> 
#>   MDES vector: 0.125, 0.125, 0.125, 0.125, 0.125
#>   nbar: 50   J: 30   K: 10   Tbar: 0.5
#>   alpha: 0.05    
#>   Level:
#>     1: R2: 0.1 (1 covariate)
#>     2: R2: 0.1 (1 covariate) ICC: 0.2    omega: 0
#>     3:   fixed effects   ICC: 0.2    omega: 0.1
#>   rho = 0.5
d_m(pp)
#> [1] "d3.2_m3ff2rc"
pump_type(pp)
#> [1] "power"
is.pumpresult(pp)
#> [1] TRUE
as.data.frame(pp)
#>    MTP D1indiv D2indiv D3indiv D4indiv D5indiv indiv.mean  min1  min2  min3
#> 1 None   0.670   0.692   0.683   0.681   0.712     0.6876    NA    NA    NA
#> 2   HO   0.527   0.544   0.541   0.541   0.557     0.5420 0.784 0.651 0.519
#>    min4 complete        SE1        SE2        SE3        SE4        SE5 df1
#> 1    NA       NA 0.05043808 0.05043808 0.05043808 0.05043808 0.05043808 279
#> 2 0.423    0.358         NA         NA         NA         NA         NA  NA
dim(pp)
#> [1]  2 18
summary(pp)
#> power result: d3.2_m3ff2rc d_m with 5 outcomes
#> 
#>   MDES vector: 0.125, 0.125, 0.125, 0.125, 0.125
#>   nbar: 50   J: 30   K: 10   Tbar: 0.5
#>   alpha: 0.05    
#>   Level:
#>     1: R2: 0.1 (1 covariate)
#>     2: R2: 0.1 (1 covariate) ICC: 0.2    omega: 0
#>     3:   fixed effects   ICC: 0.2    omega: 0.1
#>   rho = 0.5
#>   MTP    D1indiv    D2indiv    D3indiv    D4indiv    D5indiv indiv.mean  min1
#>  None      0.670      0.692      0.683      0.681      0.712      0.688      
#>    SE ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )                  
#>    HO      0.527      0.544      0.541      0.541      0.557      0.542 0.784
#>   min2  min3  min4 complete df1
#>                             279
#>                                
#>  0.651 0.519 0.423    0.358    
#>  0.005 <= MCSE <= 0.008
#>  (tnum = 1000)
transpose_power_table(pp)
#> power result: d3.2_m3ff2rc d_m with 5 outcomes
#>                 power  None    HO
#>  individual outcome 1 0.670 0.527
#>                    SE (  )       
#>  individual outcome 2 0.692 0.544
#>  individual outcome 3 0.683 0.541
#>  individual outcome 4 0.681 0.541
#>  individual outcome 5 0.712 0.557
#>       mean individual 0.688 0.542
#>             1-minimum       0.784
#>             2-minimum       0.651
#>             3-minimum       0.519
#>             4-minimum       0.423
#>              complete       0.358
#>  0.005 <= MCSE <= 0.008

J <- pump_sample(d_m = "d2.1_m2fc",
  MTP = 'HO', power.definition = 'D1indiv',
  typesample = 'J', target.power = 0.7,
  nbar = 50, M = 3, MDES = 0.125,
  Tbar = 0.5, alpha = 0.05, numCovar.1 = 1,
  R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2, tnum = 1000)
  
search_path(J)
#>    step MTP target.power       pt         dx    w     power       delta
#> 1     0  HO          0.7 27.00000         NA  100 0.5700000 -0.13000000
#> 2     0  HO          0.7 29.35261         NA  100 0.7200000  0.02000000
#> 3     0  HO          0.7 31.80348         NA  100 0.6600000 -0.04000000
#> 4     0  HO          0.7 34.35261         NA  100 0.8100000  0.11000000
#> 5     0  HO          0.7 37.00000         NA  100 0.6900000 -0.01000000
#> 6     1  HO          0.7 27.79899 0.05473879  110 0.5363636 -0.16363636
#> 7     2  HO          0.7 28.81738 0.32866631  121 0.6363636 -0.06363636
#> 8     3  HO          0.7 28.95949 0.24735533  133 0.6691729 -0.03082707
#> 9     4  HO          0.7 29.12178 0.12084406  146 0.6506849 -0.04931507
#> 10    5  HO          0.7 29.34356 0.07113758  161 0.6521739 -0.04782609
#> 11    6  HO          0.7 29.59865 0.04465419  177 0.6384181 -0.06158192
#> 12    7  HO          0.7 30.66521 0.02110859 1000 0.6760000 -0.02400000
#> 13    8  HO          0.7 32.52104 0.01391496  215 0.7767442  0.07674419
#> 14    9  HO          0.7 30.90869 0.02177979 1000 0.7000000  0.00000000
#> 15    9  HO          0.7 30.90869 0.02177979 4000 0.6735000 -0.02650000
#> 16   10  HO          0.7 31.81268 0.01698794  261 0.7279693  0.02796935
#> 17   11  HO          0.7 31.73526 0.01474169  287 0.6759582 -0.02404181
#> 18   12  HO          0.7 31.87035 0.01417985 1000 0.7140000  0.01400000
#> 19   13  HO          0.7 31.73912 0.01478548  348 0.6810345 -0.01896552
#> 20   14  HO          0.7 31.83207 0.01439764 1000 0.7100000  0.01000000
#> 21   15  HO          0.7 31.75612 0.01474307  421 0.7410926  0.04109264
#> 22   16  HO          0.7 31.58274 0.01553799 1000 0.6990000 -0.00100000
#> 23   16  HO          0.7 31.58274 0.01553799 4000 0.6985000 -0.00150000
power_curve(J)   
#>   step       pt    w MTP target.power  power
#> 1    0 27.00000 2000  HO          0.7 0.6185
#> 2    0 30.94638 2000  HO          0.7 0.6765
#> 3    0 35.16185 2000  HO          0.7 0.7605
#> 4    0 39.64638 2000  HO          0.7 0.7955
#> 5    0 44.40000 2000  HO          0.7 0.8720
```
