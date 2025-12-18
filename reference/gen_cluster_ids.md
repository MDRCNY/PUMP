# Generates school and district assignments (simulation function)

Generates simple default schools and districts IDs for individual
students for the purpose of simulations. This assumes equal sized
schools in equal sized districts.

This function is beyond the main scope of calculating power, and is
instead used for simulating data. For more info on use, see the
simulation vignette.

## Usage

``` r
gen_cluster_ids(nbar, J, K)
```

## Arguments

- nbar:

  scalar; number of individuals per school.

- J:

  scalar; number of schools per district.

- K:

  scalar; number of districts.

## Value

list; school and district assignments (S.id, D.id) for each individual.
