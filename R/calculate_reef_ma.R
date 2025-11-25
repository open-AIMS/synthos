#' Calculate Reef-level Macroalgae Cover (MA)
#'
#' Calculates reef-level macroalgae cover (MA) by combining reef-level hard coral (HCC) 
#' and soft coral (SC) cover. Macroalgae is assumed to occupy the remaining available 
#' space, and the function returns reef-level MA values both as a table and an sf object.
#'
#' @title Calculate Reef-level Macroalgae Cover
#'
#' @param data_reefs_hcc A data.frame containing reef-level hard coral cover (`Value`).
#' @param data_reefs_sc A data.frame containing reef-level soft coral cover (`Value`).
#' @param data_reefs_sf An sf object of the reef sample points.
#' @param reefs_poly_sf An sf polygon object of reef boundaries.
#'
#' @return A list containing:
#'   \itemize{
#'     \item `data_reefs_ma` – Reef-level macroalgae cover in long format with `Year` and `Value`.
#'     \item `data_reefs_pts_ma_sf` – Reef-level macroalgae cover as an sf object.
#'   }
#'
#' @author Murray
#' @export
calculate_reef_ma <- function(data_reefs_hcc, data_reefs_sc, data_reefs_sf, reefs_poly_sf) {
  data_reefs_ma <- data_reefs_hcc |>
    dplyr::rename(HCC = Value) |>
    dplyr::full_join(data_reefs_sc |> dplyr::rename(SC = Value)) |>
    dplyr::mutate(
      Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
      MA = Total_Avail,
      Value = qlogis(MA)
    ) |>
    dplyr::select(-HCC, -SC, -Total_Avail, -MA)

  data_reefs_pts_ma_sf <- data_reefs_ma |>
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    sf::st_set_crs(st_crs(data_reefs_sf))
  sf::sf_use_s2(FALSE) |> suppressMessages()
  data_reefs_pts_ma_sf <- data_reefs_pts_ma_sf |>
    sf::st_intersection(reefs_poly_sf)
  sf::sf_use_s2(TRUE) |> suppressMessages()

  list(
    data_reefs_ma =  data_reefs_ma,
    data_reefs_pts_ma_sf = data_reefs_pts_ma_sf)
}
