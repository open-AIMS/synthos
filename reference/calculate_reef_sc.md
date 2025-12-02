# Calculate Reef-level Soft Coral Cover

Calculate Reef-level Soft Coral Cover (SC)

## Usage

``` r
calculate_reef_sc(
  spatial_grid,
  spde,
  all_effects_sc,
  data_reefs_df,
  data_reefs_sf,
  reefs_poly_sf,
  config_sp
)
```

## Arguments

- spatial_grid:

  An sfc POINT object representing the spatial grid.

- spde:

  A list containing the SPDE mesh, SPDE object, Q matrix, and A matrix.

- all_effects_sc:

  A matrix containing the synthetic soft coral cover field.

- data_reefs_df:

  A data.frame of reef sample points with columns `Longitude` and
  `Latitude`.

- data_reefs_sf:

  An sf object of the reef sample points.

- reefs_poly_sf:

  An sf polygon object of reef boundaries.

- config_sp:

  A list containing config_spuration parameters, including `years`.

## Value

A list containing:

- `data_reefs_sample_sc` – SC values at the sample-level.

- `data_reefs_sc` – Reef-level SC values in long format with `Year` and
  `Value`.

- `data_reefs_pts_sc_sf` – Reef-level SC values as an sf object.

## Details

Calculates reef-level soft coral cover (SC) by projecting the spatially
explicit synthetic soft coral cover field onto reef polygons. This
function aggregates cover estimates from the spatial grid to
reef-specific samples and returns the data in both tabular and sf
formats.

## Author

Murray
