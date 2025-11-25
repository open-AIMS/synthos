#' Fine-Scale Fixed Sampling Design
#'
#' Generates a fine-scale sampling hierarchy (transects, frames, and points)
#' based on a large-scale fixed sampling design. Returns an `sf` object with
#' simulated benthic cover values at the transect level, incorporating
#' site-level and transect-level random effects for HCC, SC, and MA.
#'
#' @title Generate Fine-Scale Sampling Design (Fixed)
#'
#' @param data_fixed_locs_sf An `sf` object representing the large-scale
#'   sampling design, containing:
#'   \itemize{
#'     \item `Reef` – unique reef identifier
#'     \item `Site` – unique site identifier
#'     \item `geometry` – coordinates for each site
#'   }
#'
#' @param config_fine A list containing fine-scale sampling and variance
#'   parameters, including:
#'   \itemize{
#'     \item `Number_of_transects_per_site` – number of transects per site
#'     \item `Depths` – number of depths (not used directly here)
#'     \item `Number_of_frames_per_transect` – frames per transect
#'     \item `Points_per_frame` – points per frame
#'     \item `hcc_site_sigma`, `hcc_transect_sigma`, `hcc_sigma` – random-effect SDs for HCC
#'     \item `sc_site_sigma`, `sc_transect_sigma`, `sc_sigma` – random-effect SDs for SC
#'     \item `ma_site_sigma`, `ma_transect_sigma`, `ma_sigma` – random-effect SDs for MA
#'   }
#'
#' @return An `sf` data frame representing fine-scale sampling observations,
#'   including reef, site, transect, coordinates, year, and simulated percentage
#'   cover values (`HCC`, `SC`, `MA`).
#'
#' @details
#' This function:
#' \itemize{
#'   \item Expands large-scale site locations into multiple transects
#'   \item Adds site-level and transect-level random effects for HCC, SC, MA
#'   \item Converts simulated link-scale values into percentage cover
#'   \item Returns a tidy fine-scale dataset ready for analysis or modelling
#' }
#'
#' @author Murray
#' @export
sampling_design_fine_scale_fixed <- function(data_fixed_locs_sf, config_fine) {
  set.seed(config_fine$seed)
  data_fixed_locs_obs <- data_fixed_locs_sf |>
    dplyr::bind_cols(data_fixed_locs_sf |>
                sf::st_coordinates() |>
                as.data.frame() |>
                dplyr::rename(Longitude = X, Latitude = Y)) |>
    sf::st_drop_geometry() |>
    as.data.frame() |>
    dplyr::group_by(Longitude, Latitude, Reef) |>
    tidyr::crossing(
      Transect = paste0("T",1:config_fine$Number_of_transects_per_site)) |>
    dplyr::group_by(Site, .add = TRUE) |>
    dplyr::mutate(
      SiteEffects_HCC = rnorm(1, 0, config_fine$hcc_site_sigma),
      SiteEffects_SC = rnorm(1, 0, config_fine$sc_site_sigma),
      SiteEffects_MA = rnorm(1, 0, config_fine$ma_site_sigma)
    ) |>
    dplyr::group_by(Transect, .add = TRUE) |>
    dplyr::mutate(
      TransectEffects_HCC = rnorm(1, 0, config_fine$hcc_transect_sigma),
      TransectEffects_SC = rnorm(1, 0, config_fine$sc_transect_sigma),
      TransectEffects_MA = rnorm(1, 0, config_fine$ma_transect_sigma)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      HCC1 = HCC + SiteEffects_HCC +
        TransectEffects_HCC +
        rnorm(dplyr::n(), 0, config_fine$hcc_sigma),
      HCC2 = 100*plogis(HCC1),
      SC1 = SC + SiteEffects_SC + TransectEffects_SC +
        rnorm(dplyr::n(), 0, config_fine$sc_sigma),
      SC2 = 100*plogis(SC1),
      MA1 = MA + SiteEffects_MA + TransectEffects_MA
      + rnorm(dplyr::n(), 0, config_fine$ma_sigma),
      MA2 = 100*plogis(MA1)
    ) |>
    dplyr::arrange(Reef, Site, Transect, Year) |>
    dplyr::select(Reef, Longitude, Latitude, Site,
      Transect, Year, HCC = HCC2, SC = SC2, MA = MA2) |>
    dplyr::mutate(Date = as.POSIXct(paste0(Year, "-01-01 14:00:00")))
  return(data_fixed_locs_obs)
}
