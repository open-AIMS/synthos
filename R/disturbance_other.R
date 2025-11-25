#' Other Disturbance Layer
#'
#' Generates a synthetic disturbance layer representing other effects (e.g., crown-of-thorns, disease)
#' by combining a temporal trend with a spatial random field projected onto a spatial grid.
#'
#' @title Other Disturbance Layer
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param config_sp A list with:
#'   \itemize{
#'     \item `years` – vector of years to simulate
#'     \item `seed` – random seed
#'   }
#'
#' @return A list with:
#'   \itemize{
#'     \item `other_effects` – matrix of spatial random field values
#'     \item `other_pts_sample` – other effects projected onto the spatial grid
#'     \item `other_pts_effects` – long-format data frame suitable for plotting
#'   }
#'
#' @author Murray
#' @export
disturbance_other <- function(spatial_grid, spde, config_sp) {
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
    sort(c("years", "seed")),
    sort(names(config_sp))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  set.seed(config_sp$seed + 1)
  other_sample <- INLA::inla.qsample(length(config_sp$years),
    spde$Q,
    seed = config_sp$seed + 1,
    constr = spde$spde$f$extraconstr
  ) |>
    suppressMessages() |>
    suppressWarnings()

  rho <- rep(0.7, length(config_sp$years))
  rho <- rbeta(length(config_sp$years), 0.2, 1)
  x <- other_sample
  for (j in 2:length(config_sp$years)) {
    x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * other_sample[, j]
  }
  other_effects <- scales::rescale(x, to = c(0, 1))
  other_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    other_effects
  )

  other_pts_effects <- other_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    ## dplyr::mutate(Year = as.numeric(Year)) # ,
    dplyr::mutate(Year = config_sp$years[as.numeric(Year)])
  ## Value=scales::rescale(Value, to=c(0,1)))

  list(
    other_effects = other_effects,
    other_pts_sample = other_pts_sample,
    other_pts_effects = other_pts_effects
  )
}
