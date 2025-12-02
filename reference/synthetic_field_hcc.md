# Synthetic Hard Coral Cover Field

Synthetic Hard Coral Cover Field

## Usage

``` r
synthetic_field_hcc(
  spatial_grid,
  all_effects_df,
  baseline_sample_hcc,
  spde,
  config_sp
)
```

## Arguments

- spatial_grid:

  An `sfc_POINT` object representing the full spatial grid.

- all_effects_df:

  A data frame containing the effects of the year.

- baseline_sample_hcc:

  A data frame containing the baseline hard coral cover sample.

- spde:

  A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
  and projection matrix `A`.

- config_sp:

  A list containing simulation config_spuration parameters.

## Value

A list with:

- `all_effects_hcc` – matrix of combined baseline and yearly effects

- `all_pts_sample_hcc` – projected synthetic field onto the spatial grid

- `all_pts_effects_hcc` – long-format data frame suitable for plotting

## Details

Creates a broad-scale synthetic field of hard coral cover by combining
the baseline field with the effects of the year. All effects are assumed
to be on the link scale (logit) and projected onto the SPDE grid.

## Author

Murray
