# Pointify polygons

Convert polygon into points

## Usage

``` r
pointify_polygons(reefs_sf)
```

## Arguments

- reefs_sf:

  A spatial object of class sf with polygons.

## Value

A list with two elements: data_reefs_sf and data_reefs_df.

## Details

Convert a polygon into points.

- rasterize the reefs frame

- convert to points (centroids of raster cells)

- filter to the values of 1

- extract coordinates

- convert to data frame

## Author

Murray

## Examples

``` r
library(sf)
library(gstat)
library(ggplot2)
config <- list(
 seed = 1,
 crs = 4326,
 model = "Exp",
 psill = 1,
 range = 15,
 nugget = 0,
 patch_threshold = 1.75,
 reef_width = 0.01
)
spatial_domain <- st_geometry(
  st_multipoint(
    x = rbind(
      c(0, -10),
      c(3, -10),
      c(10, -20),
      c(1, -21),
      c(2, -16),
      c(0, -10)
    )
  )
) |>
  st_set_crs(config$crs) |>
  st_cast("POLYGON")
set.seed(config$seed)
spatial_grid <- spatial_domain |>
  st_set_crs(NA) |>
  st_sample(size = 10000, type = "regular") |>
  st_set_crs(config$crs)
simulated_field <- generate_field(spatial_grid, config)
#> [using unconditional Gaussian simulation]
simulated_patches <- generate_patches(simulated_field, config)
#> Spherical geometry (s2) switched off
#> Spherical geometry (s2) switched on
simulated_reefs <- generate_reefs(simulated_patches, config)
#> Spherical geometry (s2) switched off
#> Spherical geometry (s2) switched on
reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
```
