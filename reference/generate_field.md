# Generate a random field

Generate a simulated field

## Usage

``` r
generate_field(spatial_grid, config_sp)
```

## Arguments

- spatial_grid:

  An sf object that represents the spatial grid

- config_sp:

  A list that contains the parameters for the variogram model. The list
  should contain the following parameters:

  - seed: an integer that sets the seed for the random number generator

  - psill: a numeric value that represents the sill of the variogram
    model

  - model: a string that represents the variogram model. The following
    models are supported: "Sph", "Exp", "Gau", "Lin", "Mat", "Ste",
    "Pen", "Hug", "Hol", "Cor", "Sphlin", "Sphexp", "Sphgaus", "Sphmat",
    "Sphste", "Sphpen", "Sphhug", "Sphhol", "Sphcor"

  - range: a numeric value that represents the range of the variogram
    model

  - nugget: a numeric value that represents the nugget of the variogram
    model

## Value

an sf object that represents the simulated field

## Details

Generate a simulated field using a variogram model and an sf object. In
this context a field represents a spatially continuous surface that
represents some sort of landscape. This could be a forest or marine
system or any other.

## Author

Murray

## Examples

``` r
library(sf)
library(gstat)
library(ggplot2)
config_sp <- list(
 seed = 1,
 crs = 4326,
 model = "Exp",
 psill = 1,
 range = 15,
 nugget = 0
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
  st_set_crs(config_sp$crs) |>
  st_cast("POLYGON")
set.seed(config_sp$seed)
spatial_grid <- spatial_domain |>
  st_set_crs(NA) |>
  st_sample(size = 10000, type = "regular") |>
  st_set_crs(config_sp$crs)
simulated_field <- generate_field(spatial_grid, config_sp)
#> [using unconditional Gaussian simulation]
simulated_field |> ggplot() + geom_sf(aes(colour = sim1))
```
