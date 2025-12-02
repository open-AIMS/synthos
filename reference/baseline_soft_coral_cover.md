# Baseline Soft Coral Cover

Baseline Soft Coral Cover

## Usage

``` r
baseline_soft_coral_cover(spatial_grid, spde, config_sp)
```

## Arguments

- spatial_grid:

  An `sfc_POINT` object representing the full spatial grid.

- spde:

  A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
  and projection matrix `A`.

- config_sp:

  A list containing config_spuration parameters including `years`.

- cover_range:

  Numeric vector of length 2 defining the broad-scale range of soft
  coral cover on the link scale. Values must be \>0 and \<1. Default
  `c(0.01, 0.1)`.

## Value

A list with:

- `baseline_sample_sc` – baseline soft coral cover sample

- `baseline_effects_sc` – matrix of baseline effects

- `baseline_pts_sample_sc` – projected baseline onto the spatial grid

- `baseline_pts_effects_sc` – long-format data frame suitable for
  plotting

## Details

Calculates the baseline spatial pattern of soft coral cover prior to
sampling. The pattern is defined as a simple sine wave (applied to
centered latitudes) and projected onto the SPDE grid. Values are on the
link (logit) scale.

## Author

Murray
