#' Synthetic Soft Coral Cover Field
#'
#' Generates a broad-scale synthetic field of soft coral cover by combining
#' baseline soft coral cover with annual effects. The resulting values are
#' on the link (logit) scale and projected onto the SPDE grid.
#'
#' @title Synthetic Soft Coral Cover Field
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param all_effects_df A data.frame containing the annual effects of disturbances.
#' @param baseline_sample_sc A matrix or data.frame containing the baseline soft coral cover.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param config A list containing configuration parameters, including `years`.
#'
#' @return A list with:
#'   \itemize{
#'     \item `all_effects_sc` – matrix of combined baseline and effects on the SPDE mesh
#'     \item `all_pts_sample_sc` – projected synthetic field onto the spatial grid
#'     \item `all_pts_effects_sc` – long-format data frame suitable for plotting
#'   }
#'
#' @author Murray
#' @export

synthetic_field_sc <- function(spatial_grid, all_effects_df, baseline_sample_sc, spde, config) {
  testthat::expect(
    inherits(all_effects_df, c("data.frame")),
    "all_effects_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(all_effects_df),
    c("Longitude", "Latitude", "Year", "Y", "Growth_SC", "Y_SC")
  )
  testthat::expect(
    inherits(baseline_sample_sc, c("data.frame")),
    "baseline_sample_sc must be a data.frame object"
  )
  testthat::expect_contains(
    names(baseline_sample_sc),
    c("Longitude", "Latitude", "Y")
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  ## Do all this on the link scale so that can use cumsum
  all_effects_sc <- all_effects_df |>
    dplyr::full_join(baseline_sample_sc |>
                dplyr::select(Longitude, Latitude, BASE_SC = Y)) |>
    dplyr::group_by(Longitude, Latitude) |>
    dplyr::mutate(SC = BASE_SC + Y_SC) |>
    dplyr::ungroup() |>
    dplyr::select(-BASE_SC, -Y_SC) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = SC
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  all_pts_sample_sc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects_sc
  )
  all_pts_effects_sc <- all_pts_sample_sc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
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
  list(all_effects_sc = all_effects_sc,
    all_pts_sample_sc = all_pts_sample_sc,
    all_pts_effects_sc = all_pts_effects_sc)
}

