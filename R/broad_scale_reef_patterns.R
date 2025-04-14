##' Calculate reef-level Hard Coral Cover
##'
##' Calculate reef-level Hard Coral Cover (HCC) from the spatial grid and SPDE model.
##' @title Calculate reef-level Hard Coral Cover 
##' @param spatial_grid 
##' @param spde 
##' @param all_effects_hcc 
##' @param data_reefs_df 
##' @param data_reefs_sf 
##' @param reefs_poly_sf 
##' @return A list containing the following elements:
##' data_reefs_sample_hcc - The sample-level HCC values
##' data_reefs_hcc - The reef-level HCC values
##' data_reefs_pts_hcc_sf - The reef-level HCC values as an sf object
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' reefs_hcc <- calculate_reef_hcc(
##'   spatial_grid,
##'   matern_projection,
##'   field_hcc$all_effects_hcc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' sf_use_s2(FALSE)
##' reefs_hcc$data_reefs_pts_hcc_sf |>
##'   st_crop(xmin = 1.8, xmax = 2.3, ymin = -20.4, ymax = -19.9) |>
##'   mutate(Y = plogis(Value)) |>
##'   ggplot() +
##'   geom_sf(aes(color = Y)) +
##'   facet_wrap(~Year, nrow = 2) +
##'   scale_color_gradientn(colors = terrain.colors(10)) +
##'   coord_sf(
##'     crs = 4236,
##'     xlim = c(1.8, 2.3),
##'     ylim = c(-20.4, -19.9)
##'   ) +
##'   theme_bw(base_size = 12) +
##'   theme(axis.title = element_blank())
##' sf_use_s2(TRUE)
##' @export
calculate_reef_hcc <- function(spatial_grid, spde, all_effects_hcc, data_reefs_df, data_reefs_sf, reefs_poly_sf) {
  
  testthat::expect(
    inherits(all_effects_hcc, c("matrix")),
    "all_effects_hcc must be a matrix object"
  )
  testthat::expect(
    inherits(data_reefs_df, c("data.frame")),
    "data_reefs_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(data_reefs_df),
    c("Longitude", "Latitude")
  )
  testthat::expect(
    inherits(data_reefs_sf, c("sf")),
    "data_reefs_sf must be a sf object"
  )
  testthat::expect(
    inherits(reefs_poly_sf, c("sf")),
    "data_reefs_sf must be a sf object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  data_reefs_sample_hcc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(data_reefs_df[, 1:2]),
    all_effects_hcc
  )
  data_reefs_hcc <- data_reefs_sample_hcc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(data_reefs_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Year = config$years[as.numeric(Year)],
      Value = Value
    )

  data_reefs_pts_hcc_sf <- data_reefs_hcc |>
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(st_crs(data_reefs_sf))
  sf::sf_use_s2(FALSE) |> suppressMessages()
  data_reefs_pts_hcc_sf <- data_reefs_pts_hcc_sf |>
    sf::st_intersection(reefs_poly_sf)
  sf::sf_use_s2(TRUE) |> suppressMessages()
  list(data_reefs_sample_hcc = data_reefs_sample_hcc,
    data_reefs_hcc =  data_reefs_hcc,
    data_reefs_pts_hcc_sf = data_reefs_pts_hcc_sf)
}


##' Calculate reef-level Soft Coral Cover
##'
##' Calculate reef-level Soft Coral Cover (SC) from the spatial grid and SPDE model.
##' @title Calculate reef-level Soft Coral Cover 
##' @param spatial_grid 
##' @param spde 
##' @param all_effects_sc 
##' @param data_reefs_df 
##' @param data_reefs_sf 
##' @param reefs_poly_sf 
##' @return A list containing the following elements:
##' data_reefs_sample_sc - The sample-level SC values
##' data_reefs_sc - The reef-level SC values
##' data_reefs_pts_sc_sf - The reef-level SC values as an sf object
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_sc$baseline_sample_sc, matern_projection)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' reefs_sc <- calculate_reef_sc(
##'   spatial_grid,
##'   matern_projection,
##'   field_sc$all_effects_sc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' sf_use_s2(FALSE)
##' reefs_sc$data_reefs_pts_sc_sf |>
##'   st_crop(xmin = 1.8, xmax = 2.3, ymin = -20.4, ymax = -19.9) |>
##'   mutate(Y = plogis(Value)) |>
##'   ggplot() +
##'   geom_sf(aes(color = Y)) +
##'   facet_wrap(~Year, nrow = 2) +
##'   scale_color_gradientn(colors = terrain.colors(10)) +
##'   coord_sf(
##'     crs = 4236,
##'     xlim = c(1.8, 2.3),
##'     ylim = c(-20.4, -19.9)
##'   ) +
##'   theme_bw(base_size = 12) +
##'   theme(axis.title = element_blank())
##' sf_use_s2(TRUE)
##' @export
calculate_reef_sc <- function(spatial_grid, spde, all_effects_sc, data_reefs_df, data_reefs_sf, reefs_poly_sf) {
  
  testthat::expect(
    inherits(all_effects_sc, c("matrix")),
    "all_effects_sc must be a matrix object"
  )
  testthat::expect(
    inherits(data_reefs_df, c("data.frame")),
    "data_reefs_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(data_reefs_df),
    c("Longitude", "Latitude")
  )
  testthat::expect(
    inherits(data_reefs_sf, c("sf")),
    "data_reefs_sf must be a sf object"
  )
  testthat::expect(
    inherits(reefs_poly_sf, c("sf")),
    "data_reefs_sf must be a sf object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  data_reefs_sample_sc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(data_reefs_df[, 1:2]),
    all_effects_sc
  )
  data_reefs_sc <- data_reefs_sample_sc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(data_reefs_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Year = config$years[as.numeric(Year)],
      Value = Value
    )

  data_reefs_pts_sc_sf <- data_reefs_sc |>
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(st_crs(data_reefs_sf))
  sf::sf_use_s2(FALSE) |> suppressMessages()
  data_reefs_pts_sc_sf <- data_reefs_pts_sc_sf |>
    sf::st_intersection(reefs_poly_sf)
  sf_use_s2(TRUE) |> suppressMessages()
  list(data_reefs_sample_sc = data_reefs_sample_sc,
    data_reefs_sc =  data_reefs_sc,
    data_reefs_pts_sc_sf = data_reefs_pts_sc_sf)
}
##' Calculate reef-level Macroalgae Cover
##'
##' Calculate reef-level Macroalgae Cover (MA) from the spatial grid and SPDE model
##' @title Calculate reef-level Macroalgae Cover 
##' @param data_reefs_hcc 
##' @param data_reefs_sc 
##' @param data_reefs_sf 
##' @param reefs_poly_sf 
##' @return A list containing the following elements:
##' data_reefs_ma - The reef-level MA values
##' data_reefs_pts_ma_sf - The reef-level MA values as an sf object
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_sc$baseline_sample_sc, matern_projection)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' reefs_hcc <- calculate_reef_hcc(
##'   spatial_grid,
##'   matern_projection,
##'   field_hcc$all_effects_hcc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_sc <- calculate_reef_sc(
##'   spatial_grid,
##'   matern_projection,
##'   field_sc$all_effects_sc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_ma <- calculate_reef_ma(
##'   reefs_hcc$data_reefs_hcc,
##'   reefs_sc$data_reefs_sc,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' sf_use_s2(FALSE)
##' reefs_ma$data_reefs_pts_ma_sf |>
##'   st_crop(xmin = 1.8, xmax = 2.3, ymin = -20.4, ymax = -19.9) |>
##'   mutate(Y = plogis(Value)) |>
##'   ggplot() +
##'   geom_sf(aes(color = Y)) +
##'   facet_wrap(~Year, nrow = 2) +
##'   scale_color_gradientn(colors = terrain.colors(10)) +
##'   coord_sf(
##'     crs = 4236,
##'     xlim = c(1.8, 2.3),
##'     ylim = c(-20.4, -19.9)
##'   ) +
##'   theme_bw(base_size = 12) +
##'   theme(axis.title = element_blank())
##' sf_use_s2(TRUE)
##' @author Murray
##' @export
calculate_reef_ma <- function(data_reefs_hcc, data_reefs_sc, data_reefs_sf, reefs_poly_sf) {
  data_reefs_ma <- data_reefs_hcc |>
    dplyr::rename(HCC = Value) |>
    dplyr::full_join(data_reefs_sc |> dplyr::rename(SC = Value)) |>
    dplyr::mutate(
      Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
      MA = Total_Avail,
      Value = qlogis(MA)
    ) |>
    dplyr::select(-HCC, -SC, -Total_Avail, -MA)

  data_reefs_pts_ma_sf <- data_reefs_ma |>
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(st_crs(data_reefs_sf))
  sf::sf_use_s2(FALSE) |> suppressMessages()
  data_reefs_pts_ma_sf <- data_reefs_pts_ma_sf |>
    sf::st_intersection(reefs_poly_sf)
  sf::sf_use_s2(TRUE) |> suppressMessages()

  list(
    data_reefs_ma =  data_reefs_ma,
    data_reefs_pts_ma_sf = data_reefs_pts_ma_sf)
}
##' Calculate reef-level disturbance levels
##'
##' Calculate reef-level disturbance levels from the spatial grid and SPDE model
##' @param spatial_grid 
##' @param spde 
##' @param all_effects_hcc 
##' @param data_reefs_df 
##' @param data_reefs_sf 
##' @param reefs_poly_sf 
##' @return A list containing the following elements:
##' data_reefs_sample_disturb - The sample-level disturbance values
##' data_reefs_disturb - The reef-level disturbance values
##' data_reefs_pts_disturb_sf - The reef-level disturbance values as an sf object
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' reefs_cyc <- calculate_reef_disturbance(
##'   spatial_grid,
##'   matern_projection,
##'   cyc$cyc_effects,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' sf_use_s2(FALSE)
##' reefs_cyc$data_reefs_pts_disturb_sf |>
##'   st_crop(xmin = 1.8, xmax = 2.3, ymin = -20.4, ymax = -19.9) |>
##'   mutate(Y = Value) |>
##'   ggplot() +
##'   geom_sf(aes(color = Y)) +
##'   facet_wrap(~Year, nrow = 2) +
##'   scale_color_gradientn(colors = terrain.colors(10)) +
##'   coord_sf(
##'     crs = 4236,
##'     xlim = c(1.8, 2.3),
##'     ylim = c(-20.4, -19.9)
##'   ) +
##'   theme_bw(base_size = 12) +
##'   theme(axis.title = element_blank())
##' sf_use_s2(TRUE)
##' @export
calculate_reef_disturbances <- function(spatial_grid, spde, all_effects_disturb, data_reefs_df, data_reefs_sf, reefs_poly_sf) {
  
  testthat::expect(
    inherits(all_effects_disturb, c("matrix")),
    "all_effects_disturb must be a matrix object"
  )
  testthat::expect(
    inherits(data_reefs_df, c("data.frame")),
    "data_reefs_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(data_reefs_df),
    c("Longitude", "Latitude")
  )
  testthat::expect(
    inherits(data_reefs_sf, c("sf")),
    "data_reefs_sf must be a sf object"
  )
  testthat::expect(
    inherits(reefs_poly_sf, c("sf")),
    "data_reefs_sf must be a sf object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  data_reefs_sample_disturb <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(data_reefs_df[, 1:2]),
    all_effects_disturb
  )
  data_reefs_disturb <- data_reefs_sample_disturb |>
    as.matrix() |>
    as.data.frame() |>
    cbind(data_reefs_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Year = config$years[as.numeric(Year)],
      Value = Value
    )

  data_reefs_pts_disturb_sf <- data_reefs_disturb |>
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(st_crs(data_reefs_sf))
  sf::sf_use_s2(FALSE) |> suppressMessages()
  data_reefs_pts_disturb_sf <- data_reefs_pts_disturb_sf |>
    sf::st_intersection(reefs_poly_sf)
  sf::sf_use_s2(TRUE) |> suppressMessages()
  list(data_reefs_sample_disturb = data_reefs_sample_disturb,
    data_reefs_disturb =  data_reefs_disturb,
    data_reefs_pts_disturb_sf = data_reefs_pts_disturb_sf)
}

##' Combine reef-level benthos data
##'
##' Combine reef-level benthos data into a single data frame
##' @title Combine reef-level benthos data
##' @param data_reefs_pts_hcc_sf
##' @param data_reefs_pts_sc_sf
##' @param data_reefs_pts_ma_sf
##' @return A data frame containing the combined reef-level benthos data
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config)
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_sc$baseline_sample_sc, matern_projection)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' reefs_hcc <- calculate_reef_hcc(
##'   spatial_grid,
##'   matern_projection,
##'   field_hcc$all_effects_hcc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_sc <- calculate_reef_sc(
##'   spatial_grid,
##'   matern_projection,
##'   field_sc$all_effects_sc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_ma <- calculate_reef_ma(
##'   reefs_hcc$data_reefs_hcc,
##'   reefs_sc$data_reefs_sc,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' benthos_reefs_pts <- combine_reef_benthos(
##'   reefs_hcc$data_reefs_pts_hcc_sf,
##'   reefs_sc$data_reefs_pts_sc_sf,
##'   reefs_ma$data_reefs_pts_ma_sf
##' )
##' @author Murray
##' @export
combine_reef_benthos <- function(data_reefs_pts_hcc_sf, data_reefs_pts_sc_sf, data_reefs_pts_ma_sf) {
  data_reefs_pts_hcc_sf |>
    dplyr::rename(HCC = Value) |>
    dplyr::bind_cols(data_reefs_pts_sc_sf |>
      dplyr::select(SC = Value) |>
      sf::st_drop_geometry()) |>
    dplyr::bind_cols(data_reefs_pts_ma_sf |>
      dplyr::select(MA = Value) |>
      sf::st_drop_geometry())
}

##' Combine reef-level disturbance data
##'
##' Combine reef-level disturbance data into the benthos data frame
##' @title Combine reef-level disturbance data with benthos data
##' @param benthos_reef_pts 
##' @param data_reefs_pts_cyc_sf
##' @param data_reefs_pts_dhw_sf
##' @param data_reefs_pts_other_sf
##' @return A data frame containing the combined reef-level benthos data
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config)
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_sc$baseline_sample_sc, matern_projection)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' reefs_hcc <- calculate_reef_hcc(
##'   spatial_grid,
##'   matern_projection,
##'   field_hcc$all_effects_hcc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_sc <- calculate_reef_sc(
##'   spatial_grid,
##'   matern_projection,
##'   field_sc$all_effects_sc,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_ma <- calculate_reef_ma(
##'   reefs_hcc$data_reefs_hcc,
##'   reefs_sc$data_reefs_sc,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_cyc <- calculate_reef_disturbances(
##'   spatial_grid,
##'   matern_projection,
##'   cyc$cyc_effects,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_dhw <- calculate_reef_disturbances(
##'   spatial_grid,
##'   matern_projection,
##'   dhw$dhw_effects,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' reefs_other <- calculate_reef_disturbances(
##'   spatial_grid,
##'   matern_projection,
##'   other$other_effects,
##'   reefs$data_reefs_df,
##'   reefs$data_reefs_sf,
##'   simulated_reefs$simulated_reefs_poly_sf
##' )
##' benthos_reefs_pts <- combine_reef_benthos(
##'   reefs_hcc$data_reefs_pts_hcc_sf,
##'   reefs_sc$data_reefs_pts_sc_sf,
##'   reefs_ma$data_reefs_pts_ma_sf
##' )
##' benthos_reefs_pts <- combine_reef_disturbances(
##'   benthos_reefs_pts,
##'   reefs_hcc$data_reefs_pts_cyc_sf,
##'   reefs_sc$data_reefs_pts_dhw_sf,
##'   reefs_ma$data_reefs_pts_other_sf
##' )
##' @author Murray
##' @export
combine_reef_disturbances <- function(benthos_reefs_pts, data_reefs_pts_cyc_sf, data_reefs_pts_dhw_sf, data_reefs_pts_other_sf) {
  benthos_reefs_pts |>
    dplyr::bind_cols(data_reefs_pts_cyc_sf |>
                       dplyr::select(CYC = Value) |>
                       sf::st_drop_geometry()) |>
    dplyr::bind_cols(data_reefs_pts_dhw_sf |>
                       dplyr::select(DHW = Value) |>
                       sf::st_drop_geometry()) |>
    dplyr::bind_cols(data_reefs_pts_other_sf |>
                       dplyr::select(OTHER = Value) |>
                       sf::st_drop_geometry())
}

##' Create synthetic reef landscape
##'
##' Create a synthetic reef landscape from the spatial grid and a set of
##' user configuration parameters.  Essentially, this is a convenience
##' wrapper around many lower level functions.
##' @title Create synthetic reef landscape
##' @param spatial_grid
##' An sf object that represents the spatial grid
##' @param config
##' A list containing the following:
##' - seed: an integer that sets the seed for the random number generator
##' - psill: a numeric value that represents the sill of the variogram model
##' - model: a string that represents the variogram model. The following models are supported: "Sph", "Exp", "Gau", "Lin", "Mat", "Ste", "Pen", "Hug", "Hol", "Cor", "Sphlin", "Sphexp", "Sphgaus", "Sphmat", "Sphste", "Sphpen", "Sphhug", "Sphhol", "Sphcor"
##' - range: a numeric value that represents the range of the variogram model
##' - nugget: a numeric value that represents the nugget of the variogram model
##' - patch_threshold: a numeric representing the threshold below which
##' - years: a numeric vector that represents the years to simulate
##' - dhw_weight: a numeric value that represents the weight of the DHW disturbance
##' - cyc_weight: a numeric value that represents the weight of the CYC disturbance
##' - other_weight: a numeric value that represents the weight of the OTHER disturbance
##' - hcc_cover_range: a numeric vector of length 2 representing the baseline range of hard coral cover (e.g. c(0.1, 0.7) equates to a range of 10% to 70% cover).
##' - hcc_growth: a numeric value that represents the growth rate of the hard coral cover
##' - sc_cover_range: a numeric vector of length 2 representing the baseline range of soft coral cover (e.g. c(0.01, 0.1) equates to a range of 0.1% to 10% cover).
##' - sc_growth: a numeric value that represents the growth rate of the soft coral cover
##' @param include_disturbances: a logical value that determines whether to include disturbances in the output file (default = TRUE)
##' @param verbose
##' A logical value that determines whether to print progress messages
##' @return benthos_reefs_pts - A data frame containing the combined reef-level benthos data
##' @author Murray
##' @examples
##' library(sf)
##' library(stars)
##' library(gstat)
##' library(INLA)
##' config <- list(
##'   seed = 1,
##'   crs = 4326,
##'   model = "Exp",
##'   psill = 1,
##'   range = 15,
##'   nugget = 0,
##'   alpha = 2,
##'   kappa = 1,
##'   variance = 1,
##'   patch_threshold = 1.75,
##'   reef_width = 0.01,
##'   years = 1:12,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_cover_range = c(0.1, 0.7),
##'   hcc_growth = 0.3,
##'   sc_cover_range = c(0.01, 0.1),
##'   sc_growth =  0.3
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' benthos_reefs_pts <- create_synthetic_reef_landscape(spatial_grid, config)
##' @export
create_synthetic_reef_landscape <- function(spatial_grid, config, include_disturbances = FALSE, verbose = FALSE) {
  
  testthat::expect(
    inherits(spatial_grid, c("sfc")),
    "spatial_grid must be an sfc object"
  )
  testthat::expect_in(
    sort(c(
      "seed",
      "model",
      "psill",
      "range",
      "nugget",
      "patch_threshold",
      "reef_width",
      "years",
      "alpha",
      "kappa",
      "variance",
      "dhw_weight",
      "cyc_weight",
      "other_weight",
      "hcc_cover_range",
      "hcc_growth",
      "sc_cover_range",
      "sc_growth"
    )),
    sort(names(config))
  )
  if (verbose) cat("Generating synthetic field\n")
  simulated_field <- generate_field(spatial_grid, config)
  if (verbose) cat("Generating synthetic patches\n")
  simulated_patches <- generate_patches(simulated_field, config)
  if (verbose) cat("Generating synthetic reefs\n")
  simulated_reefs <- generate_reefs(simulated_patches, config)
  if (verbose) cat("Generating SPDE\n")
  matern_projection <- create_spde(spatial_grid, config)
  if (verbose) cat("Generating DHW layer\n")
  dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
  if (verbose) cat("Generating cyclone layer\n")
  cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
  if (verbose) cat("Generating other disturbance layer\n")
  other <- disturbance_other(spatial_grid, matern_projection, config)
  if (verbose) cat("Combine all effect layers\n")
  all_disturbance_effects <- disturbance_all(
    spatial_grid,
    dhw_effects = dhw$dhw_effects,
    cyc_effects = cyc$cyc_effects,
    other_effects = other$other_effects,
    matern_projection,
    config) 
  if (verbose) cat("Generate baseline hard coral cover\n")
  baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection,
    cover_range = config$hcc_cover_range)
  if (verbose) cat("Generate synthetic hard coral cover\n")
  field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
    baseline_hcc$baseline_sample_hcc, matern_projection)
  if (verbose) cat("Generate baseline soft coral cover\n")
  baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection,
    cover_range = config$sc_cover_range)
  if (verbose) cat("Generate synthetic soft coral cover\n")
  field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
    baseline_sc$baseline_sample_sc, matern_projection)
  if (verbose) cat("Pointify polygons\n")
  reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
  if (verbose) cat("Calculate reef hard coral cover\n")
  reefs_hcc <- calculate_reef_hcc(
    spatial_grid,
    matern_projection,
    field_hcc$all_effects_hcc,
    reefs$data_reefs_df,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf
  )
  if (verbose) cat("Calculate reef soft coral cover\n")
  reefs_sc <- calculate_reef_sc(
    spatial_grid,
    matern_projection,
    field_sc$all_effects_sc,
    reefs$data_reefs_df,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf
  )
  if (verbose) cat("Calculate reef macroalgae cover\n")
  reefs_ma <- calculate_reef_ma(
    reefs_hcc$data_reefs_hcc,
    reefs_sc$data_reefs_sc,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf
  )
  if (verbose) cat("Combine reef-level benthos data\n")
  benthos_reefs_pts <- combine_reef_benthos(
    reefs_hcc$data_reefs_pts_hcc_sf,
    reefs_sc$data_reefs_pts_sc_sf,
    reefs_ma$data_reefs_pts_ma_sf
  )
  if (include_disturbances) {
    if (verbose) cat("Calculate reef disturbances\n")
    reefs_cyc <- calculate_reef_disturbances(
      spatial_grid,
      matern_projection,
      cyc$cyc_effects,
      reefs$data_reefs_df,
      reefs$data_reefs_sf,
      simulated_reefs$simulated_reefs_poly_sf
    )
    reefs_dhw <- calculate_reef_disturbances(
      spatial_grid,
      matern_projection,
      dhw$dhw_effects,
      reefs$data_reefs_df,
      reefs$data_reefs_sf,
      simulated_reefs$simulated_reefs_poly_sf
    )
    reefs_other <- calculate_reef_disturbances(
      spatial_grid,
      matern_projection,
      other$other_effects,
      reefs$data_reefs_df,
      reefs$data_reefs_sf,
      simulated_reefs$simulated_reefs_poly_sf
    )
    benthos_reefs_pts <- combine_reef_disturbances(
      benthos_reefs_pts,
      reefs_cyc$data_reefs_pts_disturb_sf,
      reefs_dhw$data_reefs_pts_disturb_sf,
      reefs_other$data_reefs_pts_disturb_sf
    )
  }
  return(benthos_reefs_pts)
}


