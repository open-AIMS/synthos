## Generating synthethic data
rm(list = ls())


library(sf)
library(stars)
library(gstat)
library(INLA)
library(ggplot2)
detach("package:synthos", unload = TRUE, character.only = TRUE)
remotes::install_github("open-AIMS/synthos@julie", force = TRUE, dependencies = FALSE)
library(synthos)

##-----------------------------#
## 1. Generate settings
##-----------------------------#

##### Generate settings
surveys <-  "random" # or  "fixed"
data_type <- "points" # or "cover"

synthos::generateSettings(nreefs = 25, nsites = 3, nyears = 15)

##-----------------------------#
## 2. Generate the spatio-temporal domain, disturbance effects and baselines
##-----------------------------#

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

##-----------------------------#
## 3. Generate sampling design
##-----------------------------#

if (surveys == "fixed") {
  locs_sf <- synthos::sampling_design_large_scale_fixed(
    benthos_reefs_pts, config_lrge
  )
  obs <- synthos::sampling_design_fine_scale_fixed(
    locs_sf, config_fine
  )
  
} else if (surveys == "random") {
  locs_sf <- synthos::sampling_design_large_scale_random(
    benthos_reefs_pts, config_lrge
  )
  obs <- synthos::sampling_design_fine_scale_random(
    locs_sf, config_fine
  )
} else {
  stop("surveys must be 'fixed' or 'random'.")
}

##-----------------------------#
## 4. Generate export table
##-----------------------------#

if (data_type == "points") {
  pts <- synthos::sampling_design_fine_scale_points(obs, config_pt)
  synthos_data <- synthos::prepare_table(pts)  
  
} else if (data_type == "cover") {
  cov <- synthos::sampling_design_fine_scale_cover(obs, config_pt)
  synthos_data <- synthos::prepare_table(cov) 
  
} else {
  stop("data_type must be 'points' or 'cover'.")
}

##-----------------------------#
## 5. Vizualisations
##-----------------------------#

# 3.1 Long-term trajectories at site level

plots <- synthos::plot_synthos(synthos_data, type = "trajectories")

purrr::walk(seq_along(plots), ~ ggsave(filename = paste0("figures/figure1.", .x, ".png"),
                                plot = plots[[.x]], width = 6, height = 10, dpi = 300))


# 3.2 Heatmaps 
plots <- synthos::plot_synthos(synthos_data, type = "heatmaps")

purrr::walk(seq_along(plots), ~ ggsave(filename = paste0("figures/figure2.", .x, ".png"),
                                plot = plots[[.x]], width = 6, height = 8, dpi = 300))

