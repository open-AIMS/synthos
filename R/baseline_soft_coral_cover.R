#' Baseline Soft Coral Cover
#'
#' Calculates the baseline spatial pattern of soft coral cover prior to sampling.
#' The pattern is defined as a simple sine wave (applied to centered latitudes)
#' and projected onto the SPDE grid. Values are on the link (logit) scale.
#'
#' @title Baseline Soft Coral Cover
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param cover_range Numeric vector of length 2 defining the broad-scale range
#'   of soft coral cover on the link scale. Values must be >0 and <1. Default `c(0.01, 0.1)`.
#' @param config_sp A list containing config_spuration parameters including `years`.
#'
#' @return A list with:
#'   \itemize{
#'     \item `baseline_sample_sc` – baseline soft coral cover sample
#'     \item `baseline_effects_sc` – matrix of baseline effects
#'     \item `baseline_pts_sample_sc` – projected baseline onto the spatial grid
#'     \item `baseline_pts_effects_sc` – long-format data frame suitable for plotting
#'   }
#'
#' @author Murray
#' @export
baseline_soft_coral_cover <- function(spatial_grid, spde, config_sp) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )

  cover_range <- config_sp$sc_cover_range
  
  testthat::expect(
    min(cover_range) > 0 & max(cover_range) < 1,
    "cover range must be between 0 and 1"
  )
  cover_range <- qlogis(cover_range)

  baseline_sample_sc <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    dplyr::mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = clong + sin(clat) + # rnorm(1,0,1) +
        1.5 * clong + -1.5 * clat
    ) |>
    dplyr::mutate(Y = scales::rescale(Y, to = cover_range))

  baseline_effects_sc <- baseline_sample_sc |>
    dplyr::select(Y) |>
    as.matrix()
  baseline_pts_sample_sc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    baseline_effects_sc
  )
  baseline_pts_effects_sc <- baseline_pts_sample_sc |>
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
    dplyr::mutate(Year = config_sp$years[as.numeric(Year)])

  list(baseline_sample_sc = baseline_sample_sc,
    baseline_effects_sc = baseline_effects_sc,
    baseline_pts_sample_sc = baseline_pts_sample_sc,
    baseline_pts_effects_sc = baseline_pts_effects_sc)
}

