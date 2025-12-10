# Create a mesh and SPDE

Create a mesh and SPDE

## Usage

``` r
create_spde(spatial_grid, config_sp)
```

## Arguments

- spatial_grid:

  A sfc POINT object representing the full spatial grid/

- config_sp:

  A list that should contains the following parameters:

  - alpha: the smoothness parameter of the Matern covariance function

  - kappa: the range parameter of the Matern covariance function

  - variance: the variance parameter of the Matern covariance function

## Value

a list containing the mesh, spde, Q, and A

## Details

This function creates a mesh and SPDE for a given spatial grid and
config_spuration. It is primarily a wrapper for the create_spde_mesh and
create_spde_matern functions.

## Author

Murray

## Examples

``` r
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(INLA)
#> Loading required package: Matrix
#> This is INLA_25.10.19 built 2025-10-19 19:10:20 UTC.
#>  - See www.r-inla.org/contact-us for how to get help.
#>  - List available models/likelihoods/etc with inla.list.models()
#>  - Use inla.doc(<NAME>) to access documentation
config_sp <- list(crs=4326, seed = 123)
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
  st_set_crs(config_sp$crs) |>
  st_cast("POLYGON")
set.seed(config_sp$seed)
spatial_grid <- spatial_domain |>
  st_set_crs(NA) |>
  st_sample(size = 10000, type = "regular") |>
  st_set_crs(config_sp$crs)
config_sp <- list(alpha = 2, kappa = 1, variance = 1)
matern_projection <- create_spde(spatial_grid, config_sp)
```
