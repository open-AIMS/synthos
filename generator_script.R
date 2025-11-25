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

# config <- list(
#   seed = 1,
#   crs = 4326,
#   model = "Exp",
#   psill = 1,
#   range = 15,
#   nugget = 0,
#   alpha = 2,
#   kappa = 1,
#   variance = 1,
#   patch_threshold = 1.75,
#   reef_width = 0.01,
#   years = 1:12,
#   dhw_weight = 0.5,
#   cyc_weight = 0.4,
#   other_weight = 0.1,
#   hcc_cover_range = c(0.1, 0.7),
#   hcc_growth = 0.3,
#   sc_cover_range = c(0.01, 0.1),
#   sc_growth =  0.3
# )

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

# config <- list(
#   years =  1:12,
#   Number_of_transects_per_site = 5,
#   Depths = 2,
#   Number_of_frames_per_transect = 100,
#   Points_per_frame = 5,
#   ## Note, the following are on the link scale
#   hcc_site_sigma = 0.5, # variability in Sites within Locations
#   hcc_transect_sigma = 0.2, # variability in Transects within Sites
#   hcc_sigma = 0.1, # random noise

#   sc_site_sigma = 0.05, # variability in Sites within Locations
#   sc_transect_sigma = 0.02, # variability in Transects within Sites
#   sc_sigma = 0.01, # random noise

#   ma_site_sigma = 0.5, # variability in Sites within Locations
#   ma_transect_sigma = 0.2, # variability in Transects within Sites
#   ma_sigma = 0.1 # random noise
# )

benthos_fixed_locs_obs <- synthos::sampling_design_fine_scale_fixed(benthos_fixed_locs_sf, config_fine)

# config <- list(
#   Depths = 2,
#   Depth_effect_multiplier = 2,
#   Number_of_transects_per_site = 5,
#   Number_of_frames_per_transect = 100,
#   Points_per_frame = 5
# )

benthos_fixed_locs_points <- synthos::sampling_design_fine_scale_points(benthos_fixed_locs_obs, config_pt)

synthetic_fixed_points <- synthos::prepare_table(benthos_fixed_locs_points)

# config <- list(
#   Depths = 2,
#   Depth_effect_multiplier = 2,
#   Number_of_transects_per_site = 5,
#   Number_of_quadrats_per_transect = 10,
#   Quad_sigma = 0.5
# )

benthos_fixed_locs_cover <- synthos::sampling_design_fine_scale_cover(benthos_fixed_locs_obs, config_pt)

synthetic_fixed_benthos_cover <- synthos::prepare_table(benthos_fixed_locs_cover)