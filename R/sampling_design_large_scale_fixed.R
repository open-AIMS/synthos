#' Large‐scale fixed sampling design
#'
#' Selects a set number of reef locations and then randomly selects a
#' fixed number of sites within each selected location.
#'
#' @title Large‐scale fixed sampling design
#'
#' @param data_reefs_pts_sf An `sf` object representing the full field.
#'   It must contain:
#'   - `Year`: numeric year
#'   - `Reef`: unique reef ID
#'   - `HCC`, `SC`, `MA`: cover values (logit scale)
#'   - `geometry`: spatial geometry
#'
#' @param config_lrge A list with:
#'   - `n_locs`: number of reefs (locations) to select
#'   - `n_sites`: number of sites per selected reef
#'   - `seed`: random seed
#'
#' @return An `sf` object representing the large-scale sampling design.
#'   Cover values remain on the logit scale.
#'
#' @author Murray
#'
#' @examples
#' config_lrge <- list(n_locs = 25, n_sites = 2, seed = 123)
#' benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(
#'   benthos_reefs_pts, config_lrge
#' )
#'
#' @export
sampling_design_large_scale_fixed <- function(data_reefs_pts_sf, config_lrge) {
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

  ## Then filter to these Reefs before selecting a single location within
  ## each of the Reefs
  data_fixed_locs_sf <- data_reefs_pts_sf |>
    dplyr::select(Reef, geometry) |>
    dplyr::distinct(.keep_all = TRUE) |>
    dplyr::group_by(Reef) |>
    dplyr::sample_n(config_lrge$n_sites) |>
    dplyr::mutate(Site = paste0("S", 1:dplyr::n())) |>
    dplyr::ungroup() |>
    sf::st_join(data_reefs_pts_sf |>
      dplyr::select(-Reef))
  return(data_fixed_locs_sf)
}


