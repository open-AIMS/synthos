#' All Disturbance Layers
#'
#' Combines all disturbance effects (DHW, cyclones, other) with coral growth to produce
#' cumulative effects per pixel. Calculations are performed on the link scale for simplicity.
#' Macroalgae are calculated as the remaining available space: `MA = Total space - HCC - SC`.
#'
#' @title All Disturbance Layers
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param dhw_effects Matrix of DHW disturbance effects.
#' @param cyc_effects Matrix of cyclone disturbance effects.
#' @param other_effects Matrix of other disturbance effects.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param config A list with:
#'   \itemize{
#'     \item `years` – vector of years to simulate
#'     \item `seed` – random seed
#'     \item `dhw_weight` – relative influence of DHW
#'     \item `cyc_weight` – relative influence of cyclones
#'     \item `other_weight` – relative influence of other disturbances
#'     \item `hcc_growth` – annual growth rate of hard coral
#'     \item `sc_growth` – annual growth rate of soft coral
#'   }
#'
#' @return A list with:
#'   \itemize{
#'     \item `disturb_effects` – combined disturbance matrix
#'     \item `all_effects_df` – long-format data frame with cumulative effects
#'     \item `all_effects` – wide-format matrix of cumulative HCC effects
#'     \item `disturb_pts_sample` – projected effects onto the spatial grid
#'     \item `disturb_pts_effects` – long-format data frame for plotting
#'   }
#'
#' @author Murray
#' @export
disturbance_all <- function(spatial_grid, dhw_effects, cyc_effects, other_effects, spde, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed", "dhw_weight", "cyc_weight", "other_weight")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  disturb_effects <-
    (config$dhw_weight * dhw_effects) +
    (config$cyc_weight * cyc_effects) +
    (config$other_weight * other_effects) |>
    as.data.frame()
  all_effects_df <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::rename(Longitude = V1, Latitude = V2) |>
    cbind(disturb_effects) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = "Year",
      names_pattern = "sample:(.*)",
      values_to = "Y"
    ) |>
    dplyr::mutate(Year = factor(Year, levels = sort(unique(as.numeric(
      as.character(Year)
    ))))) |>
    dplyr::group_by(Longitude, Latitude) |>
    dplyr::mutate(
      Growth_HCC = config$hcc_growth, ## Add growth onto this
      Growth_SC = config$sc_growth,
      Y_HCC = cumsum(-Y + Growth_HCC), ## cumsum on link scale will accumulate effects
      Y_SC = cumsum(-Y + Growth_SC)
    )
  all_effects <- all_effects_df |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = Y_HCC
    )

  ## Project onto the spatial grid
  disturb_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects |>
      dplyr::ungroup() |> 
      dplyr::select(-Longitude, -Latitude) |>
      as.matrix()
  )
  disturb_pts_effects <- disturb_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = config$years[as.numeric(Year)])
  list(
    disturb_effects = disturb_effects,
    all_effects_df = all_effects_df,
    all_effects =  all_effects,
    disturb_pts_sample = disturb_pts_sample,
    disturb_pts_effects = disturb_pts_effects
  )
}
