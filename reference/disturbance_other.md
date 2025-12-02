# Other Disturbance Layer

Other Disturbance Layer

## Usage

``` r
disturbance_other(spatial_grid, spde, config_sp)
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

- `other_effects` – matrix of spatial random field values

- `other_pts_sample` – other effects projected onto the spatial grid

- `other_pts_effects` – long-format data frame suitable for plotting

## Details

Generates a synthetic disturbance layer representing other effects
(e.g., crown-of-thorns, disease) by combining a temporal trend with a
spatial random field projected onto a spatial grid.

## Author

Murray
