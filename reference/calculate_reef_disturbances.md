# Calculate Reef-level Disturbance Levels

Calculate Reef-level Disturbance Levels

## Usage

``` r
calculate_reef_disturbances(
  spatial_grid,
  spde,
  all_effects_disturb,
  data_reefs_df,
  data_reefs_sf,
  reefs_poly_sf,
  config_sp
)
```

## Arguments

- spatial_grid:

  An sf object of spatial points representing the full grid.

- spde:

  A list containing the SPDE mesh and related projection matrices.

- all_effects_disturb:

  A matrix of disturbance effects on the link scale.

- data_reefs_df:

  A data.frame of reef sample points (`Longitude`, `Latitude`).

- data_reefs_sf:

  An sf object of the reef sample points.

- reefs_poly_sf:

  An sf polygon object of reef boundaries.

- config_sp:

  A list containing configuration options, including `years`.

## Value

A list containing:

- `data_reefs_sample_disturb` – Sample-level disturbance values.

- `data_reefs_disturb` – Reef-level disturbance values in long format
  with `Year` and `Value`.

- `data_reefs_pts_disturb_sf` – Reef-level disturbance values as an sf
  object.

## Details

Calculates reef-level disturbance levels by projecting the spatial grid
disturbance field (e.g., thermal stress, cyclones, other impacts) onto
reef sample points and aggregating values for each reef. Returns both
tabular and sf representations.

## Author

Murray
