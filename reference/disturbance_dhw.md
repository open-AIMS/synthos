# Degree Heating Weeks Disturbance Layer

Degree Heating Weeks (DHW) Disturbance Layer

## Usage

``` r
disturbance_dhw(spatial_grid, spde, config_sp)
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

- `dhw_temporal` – data frame of the DHW temporal trend

- `dhw_effects` – matrix of spatial random field values

- `dhw_pts_sample` – DHW effects projected onto the grid

- `dhw_pts_effects_df` – long-format data frame for plotting

## Details

Generates a synthetic DHW disturbance layer by combining a temporal
trend with a spatial random field projected onto a spatial grid.

## Author

Murray
