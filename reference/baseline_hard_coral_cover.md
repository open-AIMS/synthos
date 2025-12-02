# Baseline Hard Coral Cover

Baseline Hard Coral Cover

## Usage

``` r
baseline_hard_coral_cover(spatial_grid, spde, config_sp)
```

## Arguments

- spatial_grid:

  An `sfc_POINT` object representing the full spatial grid.

- spde:

  A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
  and projection matrix `A`.

- cover_range:

  Numeric vector of length 2 defining the range of coral cover (between
  0 and 1, excluding 0 and 1). Default is `c(0.1, 0.7)` representing 10%
  to 70% cover.

## Value

A list with:

- `baseline_sample_hcc` – data frame of baseline values on the SPDE mesh

- `baseline_effects_hcc` – matrix of baseline effects

- `baseline_pts_sample_hcc` – baseline effects projected onto the
  spatial grid

- `baseline_pts_effects_hcc` – long-format data frame suitable for
  plotting

## Details

Calculates baseline hard coral cover for the year prior to sampling. The
spatial pattern is defined as a simple sine wave applied to centered
latitudes, optionally rotated, and projected onto the SPDE grid. Values
are on the expected link scale (logit) and rescaled to the nominated
cover range.

## Author

Murray
