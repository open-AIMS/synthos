#' Synthetic Hard Coral Cover Field
#'
#' Creates a broad-scale synthetic field of hard coral cover by combining
#' the baseline field with the effects of the year. All effects are assumed
#' to be on the link scale (logit) and projected onto the SPDE grid.
#'
#' @title Synthetic Hard Coral Cover Field
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param all_effects_df A data frame containing the effects of the year.
#' @param baseline_sample_hcc A data frame containing the baseline hard coral cover sample.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param config A list containing simulation configuration parameters.
#'
#' @return A list with:
#'   \itemize{
#'     \item `all_effects_hcc` – matrix of combined baseline and yearly effects
#'     \item `all_pts_sample_hcc` – projected synthetic field onto the spatial grid
#'     \item `all_pts_effects_hcc` – long-format data frame suitable for plotting
#'   }
#'
#' @author Murray
#' @export
#' 
synthetic_field_hcc <- function(spatial_grid, all_effects_df, baseline_sample_hcc, spde, config) {
  testthat::expect(
    inherits(all_effects_df, c("data.frame")),
    "all_effects_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(all_effects_df),
    c("Longitude", "Latitude", "Year", "Y", "Growth_HCC", "Y_HCC")
  )
  testthat::expect(
    inherits(baseline_sample_hcc, c("data.frame")),
    "baseline_sample_hcc must be a data.frame object"
  )
  testthat::expect_contains(
    names(baseline_sample_hcc),
    c("Longitude", "Latitude", "Y")
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  ## Do all this on the link scale so that can use cumsum
  all_effects_hcc <- all_effects_df |>
    dplyr::full_join(baseline_sample_hcc |>
                dplyr::select(Longitude, Latitude, BASE_HCC = Y)) |>
    dplyr::group_by(Longitude, Latitude) |>
    dplyr::mutate(HCC = BASE_HCC + Y_HCC) |>
    dplyr::ungroup() |>
    dplyr::select(-BASE_HCC, -Y_HCC) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = HCC
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  all_pts_sample_hcc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects_hcc
  )
  all_pts_effects_hcc <- all_pts_sample_hcc |>
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
  list(all_effects_hcc = all_effects_hcc,
    all_pts_sample_hcc = all_pts_sample_hcc,
    all_pts_effects_hcc = all_pts_effects_hcc)
}
