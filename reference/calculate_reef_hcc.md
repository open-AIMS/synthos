# Calculate Reef-level Hard Coral Cover

Calculate Reef-level Hard Coral Cover (HCC)

## Usage

``` r
calculate_reef_hcc(
  spatial_grid,
  spde,
  all_effects_hcc,
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

- all_effects_hcc:

  A matrix containing the synthetic hard coral cover field.

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

- `data_reefs_sample_hcc` – HCC values at the sample-level.

- `data_reefs_hcc` – Reef-level HCC values in long format with `Year`
  and `Value`.

- `data_reefs_pts_hcc_sf` – Reef-level HCC values as an sf object.

## Details

Calculates reef-level hard coral cover (HCC) by projecting the spatially
explicit synthetic hard coral cover field onto reef polygons. This
function aggregates cover estimates from the spatial grid to
reef-specific samples and returns the data in both tabular and sf
formats.

## Author

Murray
