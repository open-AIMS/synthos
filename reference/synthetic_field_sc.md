# Synthetic Soft Coral Cover Field

Synthetic Soft Coral Cover Field

## Usage

``` r
synthetic_field_sc(
  spatial_grid,
  all_effects_df,
  baseline_sample_sc,
  spde,
  config_sp
)
```

## Arguments

- spatial_grid:

  An `sfc_POINT` object representing the full spatial grid.

- all_effects_df:

  A data.frame containing the annual effects of disturbances.

- baseline_sample_sc:

  A matrix or data.frame containing the baseline soft coral cover.

- spde:

  A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
  and projection matrix `A`.

- config_sp_sp_sp:

  A list containing config_sp_sp_spuration parameters, including
  `years`.

## Value

A list with:

- `all_effects_sc` – matrix of combined baseline and effects on the SPDE
  mesh

- `all_pts_sample_sc` – projected synthetic field onto the spatial grid

- `all_pts_effects_sc` – long-format data frame suitable for plotting

## Details

Generates a broad-scale synthetic field of soft coral cover by combining
baseline soft coral cover with annual effects. The resulting values are
on the link (logit) scale and projected onto the SPDE grid.

## Author

Murray
