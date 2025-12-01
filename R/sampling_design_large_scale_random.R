#' Large‐scale Randomised Sampling Design
#'
#' Randomly selects reef locations and a fixed number of sites within each
#' selected location. Valid survey years are then sampled per reef using
#' temporal constraints (via `sample_years_with_condition()`).
#'
#' @title Large‐scale Randomised Sampling Design
#'
#' @param data_reefs_pts_sf An `sf` object representing the full field. Must contain:
#' \itemize{
#'   \item `Year` — survey year (numeric)
#'   \item `Reef` — unique reef identifier
#'   \item `HCC`, `SC`, `MA` — benthic cover values (logit scale)
#'   \item `geometry` — spatial geometry
#' }
#'
#' @param config_lrge A list with:
#' \itemize{
#'   \item `n_locs` — number of reef locations to select
#'   \item `n_sites` — number of sites per selected reef
#'   \item `seed` — random seed for reproducibility
#' }
#'
#' @return An `sf` object representing the large-scale randomised sampling design.
#' Includes selected reefs, sampled sites, valid selected years, and all
#' associated cover values (still on the logit scale).
#'
#' @details
#' The function:
#' \itemize{
#'   \item Randomly samples reef locations
#'   \item Randomly selects a fixed number of sites within each reef
#'   \item Uses `sample_years_with_condition()` to select years per reef
#'   \item Filters the dataset to retain only valid reef–year combinations
#' }
#'
#' @author Julie
#' @export
sampling_design_large_scale_random <- function(data_reefs_pts_sf, config_lrge) {
  testthat::expect(
    inherits(data_reefs_pts_sf, c("sf")),
    "data_reefs_pts_sf must be an sf object"
  )
  testthat::expect_in(
    sort(c(
      "n_sites"
    )),
    sort(names(config_lrge))
  )
 # set.seed(config_lrge$seed)

  ## Then filter to these Reefs before selecting a single location within
  ## each of the Reefs
  data_random_locs_sf <- data_reefs_pts_sf |>
    dplyr::select(Reef, geometry) |>
    dplyr::distinct(.keep_all = TRUE) |>
    dplyr::group_by(Reef) |>
    dplyr::sample_n(config_lrge$n_sites) |>
    dplyr::mutate(Site = paste0("S", 1:dplyr::n())) |>
    dplyr::ungroup() |>
    sf::st_join(data_reefs_pts_sf |>
      dplyr::select(-Reef))

# Sample valid years per reef
reef_years <- data_random_locs_sf |>
  dplyr::group_by(Reef) |>
  dplyr::summarise(selected_years = list(synthos::sample_years_with_condition(Year)), .groups = "drop") |>
  tidyr::unnest(cols = c(selected_years)) |>
  dplyr::mutate(reef_year = paste(Reef, selected_years, sep=""))

data_random_locs_sf <- data_random_locs_sf |>
  dplyr::mutate(reef_year = paste(Reef, Year, sep="")) |>
  dplyr::filter(reef_year %in% reef_years$reef_year) |>
  dplyr::select(!reef_year)

  return(data_random_locs_sf)
}
