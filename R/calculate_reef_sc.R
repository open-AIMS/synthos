#' Calculate Reef-level Soft Coral Cover (SC)
#'
#' Calculates reef-level soft coral cover (SC) by projecting the spatially explicit 
#' synthetic soft coral cover field onto reef polygons. This function aggregates 
#' cover estimates from the spatial grid to reef-specific samples and returns the 
#' data in both tabular and sf formats.
#'
#' @title Calculate Reef-level Soft Coral Cover
#'
#' @param spatial_grid An sfc POINT object representing the spatial grid.
#' @param spde A list containing the SPDE mesh, SPDE object, Q matrix, and A matrix.
#' @param all_effects_sc A matrix containing the synthetic soft coral cover field.
#' @param data_reefs_df A data.frame of reef sample points with columns `Longitude` and `Latitude`.
#' @param data_reefs_sf An sf object of the reef sample points.
#' @param reefs_poly_sf An sf polygon object of reef boundaries.
#' @param config_sp A list containing config_spuration parameters, including `years`.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `data_reefs_sample_sc` – SC values at the sample-level.
#'     \item `data_reefs_sc` – Reef-level SC values in long format with `Year` and `Value`.
#'     \item `data_reefs_pts_sc_sf` – Reef-level SC values as an sf object.
#'   }
#'
#' @author Murray
#' @export
calculate_reef_sc <- function(spatial_grid, spde, all_effects_sc, data_reefs_df, data_reefs_sf, reefs_poly_sf, config_sp) {
  
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
      Year = config_sp$years[as.numeric(Year)],
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
