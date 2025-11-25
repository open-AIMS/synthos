#' Calculate Reef-level Disturbance Levels
#'
#' Calculates reef-level disturbance levels by projecting the spatial grid disturbance field
#' (e.g., thermal stress, cyclones, other impacts) onto reef sample points and aggregating
#' values for each reef. Returns both tabular and sf representations.
#'
#' @title Calculate Reef-level Disturbance Levels
#'
#' @param spatial_grid An sf object of spatial points representing the full grid.
#' @param spde A list containing the SPDE mesh and related projection matrices.
#' @param all_effects_disturb A matrix of disturbance effects on the link scale.
#' @param data_reefs_df A data.frame of reef sample points (`Longitude`, `Latitude`).
#' @param data_reefs_sf An sf object of the reef sample points.
#' @param reefs_poly_sf An sf polygon object of reef boundaries.
#' @param config_sp A list containing configuration options, including `years`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `data_reefs_sample_disturb` – Sample-level disturbance values.
#'     \item `data_reefs_disturb` – Reef-level disturbance values in long format with `Year` and `Value`.
#'     \item `data_reefs_pts_disturb_sf` – Reef-level disturbance values as an sf object.
#'   }
#'
#' @author Murray
#' @export
calculate_reef_disturbances <- function(spatial_grid, spde, all_effects_disturb, data_reefs_df, data_reefs_sf, reefs_poly_sf, config_sp) {
  
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
      Year = config_sp$years[as.numeric(Year)],
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
