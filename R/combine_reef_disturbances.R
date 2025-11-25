#' Combine Reef-level Disturbance Data
#'
#' Combines reef-level disturbance data (CYC, DHW, OTHER) with the benthos data frame,
#' preserving reef coordinates and removing geometry columns from disturbance data.
#'
#' @title Combine Reef-level Disturbance Data with Benthos
#'
#' @param benthos_reefs_pts A data frame containing reef-level benthos data (HCC, SC, MA).
#' @param data_reefs_pts_cyc_sf An sf object containing reef-level cyclonic disturbance values.
#' @param data_reefs_pts_dhw_sf An sf object containing reef-level degree heating week (DHW) values.
#' @param data_reefs_pts_other_sf An sf object containing reef-level other disturbance values.
#'
#' @return A data.frame containing reef-level benthos and disturbance values with columns:
#'         `HCC`, `SC`, `MA`, `CYC`, `DHW`, `OTHER`, along with reef coordinates.
#'
#' @author Murray
#' @export
combine_reef_disturbances <- function(benthos_reefs_pts, data_reefs_pts_cyc_sf, data_reefs_pts_dhw_sf, data_reefs_pts_other_sf) {
  benthos_reefs_pts |>
    dplyr::bind_cols(data_reefs_pts_cyc_sf |>
                       dplyr::select(CYC = Value) |>
                       sf::st_drop_geometry()) |>
    dplyr::bind_cols(data_reefs_pts_dhw_sf |>
                       dplyr::select(DHW = Value) |>
                       sf::st_drop_geometry()) |>
    dplyr::bind_cols(data_reefs_pts_other_sf |>
                       dplyr::select(OTHER = Value) |>
                       sf::st_drop_geometry())
}
