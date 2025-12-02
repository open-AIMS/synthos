# Cyclones Disturbance Layer

Cyclones Disturbance Layer

## Usage

``` r
disturbance_cyc(spatial_grid, spde, config_sp)
```

## Arguments

- spatial_grid:

  An `sfc_POINT` object representing the full spatial grid.

- spde:

  A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
  and projection matrix `A`.

- config_sp:

  A list with:

  - `years` – vector of years to simulate

  - `seed` – random seed

## Value

A list with:

- `cyc_effects` – matrix of spatial random field values

- `cyc_effects_df` – long-format data frame with `Longitude`,
  `Latitude`, `Year`, and `Value`

- `cyc_pts_sample` – cyclone effects projected onto the spatial grid

- `cyc_pts_effects` – long-format data frame suitable for plotting

## Details

Generates a synthetic cyclone disturbance layer by combining a temporal
trend with a spatial random field projected onto a spatial grid.

## Author

Murray
