#' Cyclones Disturbance Layer
#'
#' Generates a synthetic cyclone disturbance layer by combining a temporal trend
#' with a spatial random field projected onto a spatial grid.
#'
#' @title Cyclones Disturbance Layer
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
#'     \item `cyc_effects` – matrix of spatial random field values
#'     \item `cyc_effects_df` – long-format data frame with `Longitude`, `Latitude`, `Year`, and `Value`
#'     \item `cyc_pts_sample` – cyclone effects projected onto the spatial grid
#'     \item `cyc_pts_effects` – long-format data frame suitable for plotting
#'   }
#'
#' @author Murray
#' @export 
disturbance_cyc <- function(spatial_grid, spde, config) {

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

  set.seed(config$seed)
  cyc <- vector("list", length(config$years))

  yrs <- 1 + (config$years - min(config$years))
  ## for (yr in config$years) {
  for (yr in yrs) {
    ## cat(paste("Year:", yr, "\n"))
    cyc_occur <- rbinom(1, 1, prob = min(0.05 * yr^2, 0.6))
    ## cat(paste("Cyclone Occurance:", cyc_occur, "\n"))
    cyc_intensity <- rbeta(1, 2, 1) |> round(2)
    ## cat(paste("Cyclone intensity:", cyc_intensity, "\n"))
    ## cyc_spatial <- spatial_grid_pts_df  |>
    lat_offset <- runif(1, 0, 5)
    cyc_spatial <- spde$mesh$loc[, 1:2] |>
      as.data.frame() |>
      dplyr::select(Longitude = V1, Latitude = V2) |>
      dplyr::mutate(
        clong = as.vector(scale(Longitude, scale = FALSE)),
        clat = as.vector(scale(Latitude, scale = FALSE)),
        Y = lat_offset + runif(1, -1, 1) * clong + runif(1, -1, 1) *
          clat + sin(clat),
        # Y= Y - runif(1,-10,10),
        Y = abs(Y),
        Y = ifelse(Y > cyc_intensity, cyc_intensity, Y),
        Y = cyc_intensity - Y,
        Value = Y * cyc_occur
      )
    cyc[[yr]] <- cyc_spatial |>
      dplyr::mutate(Year = yr)
  }
  cyc <- do.call("rbind", cyc)
  cyc_effects_df <- cyc |>
    dplyr::mutate(Value = scales::rescale(Value, to = c(0, 1)))

  cyc_effects <- cyc_effects_df |>
    dplyr::select(-clong, -clat, -Y) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = Value
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  cyc_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    cyc_effects
  )

  cyc_pts_effects <- cyc_pts_sample |>
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
    cyc_effects = cyc_effects,
    cyc_effects_df = cyc_effects_df,
    cyc_pts_sample = cyc_pts_sample,
    cyc_pts_effects = cyc_pts_effects
  )
}
