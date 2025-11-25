#' Convert sf Geometry to a Coordinate Data Frame
#'
#' Extracts point coordinates from an sf geometry column (`sfc`) and returns
#' them as a tidy data frame with renamed longitude and latitude columns.
#'
#' @title Convert sf Geometry to Data Frame
#'
#' @param spatial_grid An sf object containing point geometries from which
#'   coordinates will be extracted.
#'
#' @return A data frame with columns `Longitude` and `Latitude`, ordered
#'   lexicographically.
#'
#' @author Murray
#' @export
spatial_grid_sfc_to_df <- function(spatial_grid) {
  spatial_grid |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    dplyr::arrange(Longitude, Latitude)
}
