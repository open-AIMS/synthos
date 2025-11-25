#' Degree Heating Weeks (DHW) Disturbance Layer
#'
#' Generates a synthetic DHW disturbance layer by combining a temporal trend
#' with a spatial random field projected onto a spatial grid.
#'
#' @title Degree Heating Weeks Disturbance Layer
#'
#' @param spatial_grid An `sfc_POINT` object representing the full spatial grid.
#' @param spde A list containing the SPDE mesh, SPDE object, precision matrix `Q`,
#'   and projection matrix `A`.
#' @param config A list with:
#'   \itemize{
#'     \item `years` – vector of years to simulate
#'     \item `seed` – random seed
#'   }
#'
#' @return A list with:
#'   \itemize{
#'     \item `dhw_temporal` – data frame of the DHW temporal trend
#'     \item `dhw_effects` – matrix of spatial random field values
#'     \item `dhw_pts_sample` – DHW effects projected onto the grid
#'     \item `dhw_pts_effects_df` – long-format data frame for plotting
#'   }
#'
#' @author Murray
#' @export
#' 
disturbance_dhw <- function(spatial_grid, spde, config) {
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
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  ## Overall temporal trend in DHW
  set.seed(config$seed)
  dhw_temporal <- data.frame(Year = config$years) |>
    dplyr::mutate(
      cYear = Year - 1, # as.vector(scale(Year, scale=FALSE)),
      Y = 0.2 * cYear + sin(cYear),
      Y = Y * rbeta(length(config$years), Y, 1),
      Y = scales::rescale(Y - min(Y), to = c(0, 5))
    )
  ## Propagate this temporal trend across a random rield with a time
  ## varying autocorrelation coefficient drawn from a beta distribution
  ## with shape parameters of 0.2 and 1
  set.seed(config$seed)
  dhw_sample <- INLA::inla.qsample(length(config$years),
    spde$Q,
    seed = config$seed,
    constr = spde$spde$f$extraconstr
  ) |>
    suppressMessages() |>
    suppressWarnings()

  rho <- rep(0.7, length(config$years))
  rho <- rbeta(length(config$years), 0.2, 1)
  x <- dhw_sample
  for (j in 2:length(config$years)) {
    x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * dhw_sample[, j]
  }
  x <- sweep(x, 2, dhw_temporal$Y, FUN = "+")
  dhw_effects <- scales::rescale(x, to = c(0, 1))
  dhw_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    dhw_effects
  )

  dhw_pts_effects_df <- dhw_pts_sample %>%
    as.matrix() %>%
    as.data.frame() %>%
    cbind(spatial_grid_pts_df) %>%
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Year = config$years[as.numeric(Year)])

  list(
    dhw_temporal = dhw_temporal,
    dhw_effects = dhw_effects,
    dhw_pts_sample = dhw_pts_sample,
    dhw_pts_effects_df = dhw_pts_effects_df
  )
}
