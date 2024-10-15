spatial_grid_sfc_to_df <- function(spatial_grid) {
  spatial_grid |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    dplyr::arrange(Longitude, Latitude)
}
