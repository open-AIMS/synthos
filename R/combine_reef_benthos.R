#' Combine Reef-level Benthos Data
#'
#' Combines reef-level benthos data (hard coral, soft coral, and macroalgae)
#' into a single data frame, preserving reef coordinates and removing geometry.
#'
#' @title Combine Reef-level Benthos Data
#'
#' @param data_reefs_pts_hcc_sf An sf object containing reef-level hard coral cover values.
#' @param data_reefs_pts_sc_sf An sf object containing reef-level soft coral cover values.
#' @param data_reefs_pts_ma_sf An sf object containing reef-level macroalgae cover values.
#'
#' @return A data.frame containing reef-level benthos values with columns `HCC`, `SC`, `MA`,
#'         along with reef coordinates (`Longitude`, `Latitude`) and any other original metadata.
#'
#' @author Murray
#' @export
combine_reef_benthos <- function(data_reefs_pts_hcc_sf, data_reefs_pts_sc_sf, data_reefs_pts_ma_sf) {
  data_reefs_pts_hcc_sf |>
    dplyr::rename(HCC = Value) |>
    dplyr::bind_cols(data_reefs_pts_sc_sf |>
      dplyr::select(SC = Value) |>
      sf::st_drop_geometry()) |>
    dplyr::bind_cols(data_reefs_pts_ma_sf |>
      dplyr::select(MA = Value) |>
      sf::st_drop_geometry())
}
