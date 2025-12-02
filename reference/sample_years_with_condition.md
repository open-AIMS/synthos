# Constrained Year Subsampling

Sample Years Under Temporal Constraints

## Usage

``` r
sample_years_with_condition(
  years,
  min_n = 2,
  max_n = 15,
  max_gap = 5,
  max_tries = 100
)
```

## Arguments

- years:

  A numeric vector of years to sample from.

- min_n:

  Minimum number of years required in the sampled subset.

- max_n:

  Maximum number of years allowed in the sampled subset.

- max_gap:

  Maximum allowable gap (in years) between consecutive selected years.

- max_tries:

  Maximum number of attempts before returning `NULL` if no valid
  combination is found.

## Value

A sorted numeric vector of sampled years satisfying the constraints, or
`NULL` if no valid combination is found after `max_tries` attempts.

## Details

Randomly selects a subset of years while enforcing constraints on
minimum sample size, maximum sample size, and the maximum allowable
temporal gap between consecutive selected years. Useful for generating
realistic temporal sampling schemes in synthetic or constrained
ecological datasets.

The function:

- Ensures sampled years fall within the specified range of sample sizes

- Enforces a temporal continuity rule (`diff(years) <= max_gap`)

- Randomly tries up to `max_tries` combinations before giving up

## Author

Julie
