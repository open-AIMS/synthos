# All Disturbance Layers

All Disturbance Layers

## Usage

``` r
disturbance_all(
  spatial_grid,
  dhw_effects,
  cyc_effects,
  other_effects,
  matern_projection,
  config_sp
)
```

## Arguments

- spatial_grid:

  An `sfc_POINT` object representing the full spatial grid.

- dhw_effects:

  Matrix of DHW disturbance effects.

- cyc_effects:

  Matrix of cyclone disturbance effects.

- other_effects:

  Matrix of other disturbance effects.

- config_sp:

  A list with:

  - `years` – vector of years to simulate

  - `seed` – random seed

  - `dhw_weight` – relative influence of DHW

  - `cyc_weight` – relative influence of cyclones

  - `other_weight` – relative influence of other disturbances

  - `hcc_growth` – annual growth rate of hard coral

  - `sc_growth` – annual growth rate of soft coral

- spde:

  A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
  and projection matrix `A`.

## Value

A list with:

- `disturb_effects` – combined disturbance matrix

- `all_effects_df` – long-format data frame with cumulative effects

- `all_effects` – wide-format matrix of cumulative HCC effects

- `disturb_pts_sample` – projected effects onto the spatial grid

- `disturb_pts_effects` – long-format data frame for plotting

## Details

Combines all disturbance effects (DHW, cyclones, other) with coral
growth to produce cumulative effects per pixel. Calculations are
performed on the link scale for simplicity. Macroalgae are calculated as
the remaining available space: `MA = Total space - HCC - SC`.

## Author

Murray
