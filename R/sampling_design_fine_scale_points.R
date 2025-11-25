#' Fine-Scale Photo-Point Sampling
#'
#' Distributes benthic cover (HCC, SC, MA) across photo-frames and point-level
#' observations. Expands fine-scale transect data to include depth effects,
#' allocates points in proportion to benthic cover, and returns a point-level
#' sampling dataset suitable for analysis or modelling.
#'
#' @title Generate Fine-Scale Photo-Point Observations
#'
#' @param data_fixed_locs_obs
#' An `sf` object representing fine-scale transect observations produced by
#' `sampling_design_fine_scale_fixed()`. Must include:
#' \itemize{
#'   \item `Reef`, `Site`, `Transect`, `Year`, `Date`
#'   \item Benthic cover values: `HCC`, `SC`, `MA`
#'   \item Transect-level coordinates
#' }
#'
#' @param config_pt A list containing fine-scale point sampling parameters:
#' \itemize{
#'   \item `Depths` – number of depths per transect
#'   \item `Depth_effect_multiplier` – magnitude of the depth effect
#'   \item `Number_of_transects_per_site` – number of transects (used for checks)
#'   \item `Number_of_frames_per_transect` – frames per transect
#'   \item `Points_per_frame` – points per frame
#'   \item `seed` – optional random seed (if provided)
#' }
#'
#' @return An `sf` object representing point-level sampling, containing:
#' \itemize{
#'   \item Reef, Site, Transect, Year, Depth, Frame, Point Number
#'   \item Group (HCC, SC, MA)
#'   \item Value – benthic cover value after depth effects
#'   \item geometry – inherited from transect-level coordinates
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Expands each transect across multiple depths
#'   \item Applies a randomised depth effect to benthic cover values
#'   \item Converts cover percentages to point allocations
#'   \item Expands into individual point observations and assigns frame/point IDs
#' }
#'
#' @author Murray
#'
#' @export
sampling_design_fine_scale_points <- function(data_fixed_locs_obs, config_pt) {
  set.seed(config_pt$seed)
  ## put on fold scale
  data_fixed_locs_obs <-
    data_fixed_locs_obs |>
    tidyr::crossing(Depth = seq(3, 10, length = config_pt$Depths)) |>
    tidyr::pivot_longer(cols = c(HCC, SC, MA),
      names_to = "Group",
      values_to = "Value") |>
    dplyr::group_by(Reef, Site, Transect, Year, Date) |>
    dplyr::mutate(Value = Value + rev(sort(config_pt$Depth_effect_multiplier *
                                      scale(rnorm(config_pt$Depths))))) |>
    dplyr::ungroup()
  ## Need to split the percentage cover into point and frames
  data_fixed_locs_obs <- data_fixed_locs_obs |>
    dplyr::group_by(Reef, Site, Transect, Year, Depth, Date) |>
    dplyr::mutate(
      Points = round(config_pt$Number_of_frames_per_transect *
                       config_pt$Points_per_frame *
                       (Value / sum(Value)), 0),
      Points = ifelse(Points < 0, 0, Points)
    ) |>
    tidyr::uncount(Points) |>
    dplyr::sample_n(dplyr::n(), replace = FALSE) |>
    dplyr::mutate(
      POINT_NO = rep_len(1:config_pt$Points_per_frame, length = dplyr::n()),
      FRAME = rep(1:config_pt$Number_of_frames_per_transect, each = config_pt$Points_per_frame, length = dplyr::n())
    ) |>
    dplyr::ungroup()
  return(data_fixed_locs_obs)
}
