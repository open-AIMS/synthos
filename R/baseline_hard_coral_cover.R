#' Baseline Hard Coral Cover
#'
#' Calculates baseline hard coral cover for the year prior to sampling.
#' The spatial pattern is defined as a simple sine wave applied to centered
#' latitudes, optionally rotated, and projected onto the SPDE grid.
#' Values are on the expected link scale (logit) and rescaled to the nominated
#' cover range.
#'
#' @title Baseline Hard Coral Cover
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param cover_range Numeric vector of length 2 defining the range of coral cover
#'   (between 0 and 1, excluding 0 and 1). Default is `c(0.1, 0.7)` representing 10% to 70% cover.
#'
#' @return A list with:
#'   \itemize{
#'     \item `baseline_sample_hcc` – data frame of baseline values on the SPDE mesh
#'     \item `baseline_effects_hcc` – matrix of baseline effects
#'     \item `baseline_pts_sample_hcc` – baseline effects projected onto the spatial grid
#'     \item `baseline_pts_effects_hcc` – long-format data frame suitable for plotting
#'   }
#'
#' @author Murray
#' @export
baseline_hard_coral_cover <- function(spatial_grid, spde, cover_range = c(0.1, 0.7)) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  testthat::expect(
    min(cover_range) > 0 & max(cover_range) < 1,
    "cover range must be between 0 and 1"
  )
  cover_range <- qlogis(cover_range)

  baseline_sample_hcc <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    dplyr::mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = clong + sin(clat) + # rnorm(1,0,1) +
        1.5 * clong + clat
    ) |>
    dplyr::mutate(Y = scales::rescale(Y, to = cover_range))

  baseline_effects_hcc <- baseline_sample_hcc |>
    dplyr::select(Y) |>
    as.matrix()
  baseline_pts_sample_hcc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    baseline_effects_hcc
  )
  baseline_pts_effects_hcc <- baseline_pts_sample_hcc |>
    cbind() |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = as.numeric(Year))

  list(baseline_sample_hcc = baseline_sample_hcc,
    baseline_effects_hcc = baseline_effects_hcc,
    baseline_pts_sample_hcc = baseline_pts_sample_hcc,
    baseline_pts_effects_hcc = baseline_pts_effects_hcc)
}
