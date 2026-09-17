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
#>  None      0.694      0.702      0.738      0.713      0.715      0.712      
#>    SE ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )                  
#>    HO      0.550      0.551      0.555      0.567      0.552      0.555 0.807
#>   min2  min3  min4 complete df1
#>                             279
#>                                
#>  0.657 0.542 0.436    0.379    
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
#> 1 None   0.694   0.702   0.738   0.713   0.715     0.7124    NA    NA    NA
#> 2   HO   0.550   0.551   0.555   0.567   0.552     0.5550 0.807 0.657 0.542
#>    min4 complete        SE1        SE2        SE3        SE4        SE5 df1
#> 1    NA       NA 0.05043808 0.05043808 0.05043808 0.05043808 0.05043808 279
#> 2 0.436    0.379         NA         NA         NA         NA         NA  NA
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
#>  None      0.694      0.702      0.738      0.713      0.715      0.712      
#>    SE ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )  ( 0.050 )                  
#>    HO      0.550      0.551      0.555      0.567      0.552      0.555 0.807
#>   min2  min3  min4 complete df1
#>                             279
#>                                
#>  0.657 0.542 0.436    0.379    
#>  0.005 <= MCSE <= 0.008
#>  (tnum = 1000)
transpose_power_table(pp)
#> power result: d3.2_m3ff2rc d_m with 5 outcomes
#>                 power  None    HO
#>  individual outcome 1 0.694 0.550
#>                    SE (  )       
#>  individual outcome 2 0.702 0.551
#>  individual outcome 3 0.738 0.555
#>  individual outcome 4 0.713 0.567
#>  individual outcome 5 0.715 0.552
#>       mean individual 0.712 0.555
#>             1-minimum       0.807
#>             2-minimum       0.657
#>             3-minimum       0.542
#>             4-minimum       0.436
#>              complete       0.379
#>  0.005 <= MCSE <= 0.008

J <- pump_sample(d_m = "d2.1_m2fc",
  MTP = 'HO', power.definition = 'D1indiv',
  typesample = 'J', target.power = 0.7,
  nbar = 50, M = 3, MDES = 0.125,
  Tbar = 0.5, alpha = 0.05, numCovar.1 = 1,
  R2.1 = 0.1, ICC.2 = 0.05, rho = 0.2, tnum = 1000)
  
search_path(J)
#>    step MTP target.power       pt         dx    w     power       delta
#> 1     0  HO          0.7 27.00000         NA  100 0.6000000 -0.10000000
#> 2     0  HO          0.7 29.35261         NA  100 0.6100000 -0.09000000
#> 3     0  HO          0.7 31.80348         NA  100 0.7200000  0.02000000
#> 4     0  HO          0.7 34.35261         NA  100 0.7600000  0.06000000
#> 5     0  HO          0.7 37.00000         NA  100 0.8500000  0.15000000
#> 6     1  HO          0.7 31.96016 0.02845146  110 0.6000000 -0.10000000
#> 7     2  HO          0.7 33.10664 0.04117286  121 0.7272727  0.02727273
#> 8     3  HO          0.7 32.84054 0.04190045  133 0.6842105 -0.01578947
#> 9     4  HO          0.7 32.93988 0.04237933  146 0.7328767  0.03287671
#> 10    5  HO          0.7 32.77112 0.04273074  161 0.6645963 -0.03540373
#> 11    6  HO          0.7 32.92222 0.04430476  177 0.7401130  0.04011299
#> 12    7  HO          0.7 32.77354 0.04493506  195 0.7589744  0.05897436
#> 13    8  HO          0.7 32.57874 0.04375803 1000 0.7100000  0.01000000
#> 14    9  HO          0.7 32.51695 0.04231615 1000 0.7160000  0.01600000
#> 15   10  HO          0.7 32.43206 0.04015805 1000 0.6900000 -0.01000000
#> 16   11  HO          0.7 32.47804 0.04176017  287 0.7317073  0.03170732
#> 17   12  HO          0.7 32.39186 0.03508235  316 0.7120253  0.01202532
#> 18   13  HO          0.7 32.36124 0.03441182  348 0.7155172  0.01551724
#> 19   14  HO          0.7 32.13187 0.02219283  383 0.7650131  0.06501305
#> 20   15  HO          0.7 32.08620 0.02766531  421 0.7173397  0.01733967
#> 21   16  HO          0.7 31.88809 0.02205892 1000 0.6890000 -0.01100000
#> 22   17  HO          0.7 32.06170 0.02737492  509 0.7308448  0.03084479
#> 23   18  HO          0.7 31.95503 0.02555570  560 0.7107143  0.01071429
#> 24   19  HO          0.7 31.91375 0.02484939 1000 0.7190000  0.01900000
#> 25   20  HO          0.7 31.81726 0.02333785 1000 0.6960000 -0.00400000
#> 26   20  HO          0.7 31.81726 0.02333785 4000 0.7050000  0.00500000
power_curve(J)   
#>   step       pt    w MTP target.power  power
#> 1    0 27.00000 2000  HO          0.7 0.6265
#> 2    0 30.94638 2000  HO          0.7 0.6915
#> 3    0 35.16185 2000  HO          0.7 0.7395
#> 4    0 39.64638 2000  HO          0.7 0.8035
#> 5    0 44.40000 2000  HO          0.7 0.8655
```
