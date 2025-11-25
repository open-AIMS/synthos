## Generating synthethic data
rm(list = ls())


library(sf)
library(stars)
library(gstat)
library(INLA)
detach("package:synthos", unload = TRUE, character.only = TRUE)
remotes::install_github("open-AIMS/synthos@julie", force = TRUE, dependencies = FALSE)
library(synthos)

##### Generate settings
synthos::generateSettings(nreefs = 25, nsites = 3, nyears = 15)

# 1. Create synthetic reef landscape
spatial_domain <- st_geometry(
  st_multipoint(
    x = rbind(
      c(0, -11),
      c(3,-11),
      c(6,-14),
      c(1,-15),
      c(2,-12),
      c(0,-11)
    )
  )
) |>
  st_set_crs(config_sp$crs) |>
  st_cast("POLYGON")

## ---- SpatialPoints
set.seed(config_sp$seed)
spatial_grid <- spatial_domain |>
  st_set_crs(NA) |>
  st_sample(size = 10000, type = "regular") |>
  st_set_crs(config_sp$crs)
sf_use_s2(FALSE)

benthos_reefs_pts <- synthos::create_synthetic_reef_landscape(spatial_grid, config_sp)
## ----end

benthos_fixed_locs_sf <- synthos::sampling_design_large_scale_fixed(benthos_reefs_pts, config_lrge)

benthos_fixed_locs_obs <- synthos::sampling_design_fine_scale_fixed(benthos_fixed_locs_sf, config_fine)

benthos_fixed_locs_points <- synthos::sampling_design_fine_scale_points(benthos_fixed_locs_obs, config_pt)

synthetic_fixed_points <- synthos::prepare_table(benthos_fixed_locs_points)

benthos_fixed_locs_cover <- synthos::sampling_design_fine_scale_cover(benthos_fixed_locs_obs, config_pt)

synthetic_fixed_benthos_cover <- synthos::prepare_table(benthos_fixed_locs_cover)
