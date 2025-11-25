
##' Fine scale photo-points
##'
##' Disperse benthos cover across the photo-frames and points 
##' @title Fine scale photo-points
##' @param data_fixed_locs_obs 
##' A sf object representing ...
##' @param config
##' A list containing:
##' - Depths: number of depths
##' - Depth_effect_multiplier: a multiplier for the depth effect
##' - Number_of_transects_per_site: number of transects per site
##' - Number_of_frames_per_transect: number of frames per transect
##' - Points_per_frame: number of points per frame
##' @return 
##' A sf object representing the fine scale sampling design with points and frames
##' @author Murray
##' @examples
##' library(sf)
##' library(stars)
##' library(gstat)
##' library(INLA)
##' config <- list(
##'   seed = 1,
##'   crs = 4326,
##'   model = "Exp",
##'   psill = 1,
##'   range = 15,
##'   nugget = 0,
##'   alpha = 2,
##'   kappa = 1,
##'   variance = 1,
##'   patch_threshold = 1.75,
##'   reef_width = 0.01,
##'   years = 1:12,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_cover_range = c(0.1, 0.7),
##'   hcc_growth = 0.3,
##'   sc_cover_range = c(0.01, 0.1),
##'   sc_growth =  0.3
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' benthos_reefs_pts <- create_synthetic_reef_landscape(spatial_grid, config)
##' config <- list(n_locs = 25, n_sites = 2, seed = 123)
##' benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(benthos_reefs_pts, config)
##' config <- list(
##'   years =  1:12,
##'   Number_of_transects_per_site = 5,
##'   Depths = 2,
##'   Number_of_frames_per_transect = 100,
##'   Points_per_frame = 5,
##'   ## Note, the following are on the link scale
##'   hcc_site_sigma = 0.5, # variability in Sites within Locations
##'   hcc_transect_sigma = 0.2, # variability in Transects within Sites
##'   hcc_sigma = 0.1, # random noise
##'   sc_site_sigma = 0.05, # variability in Sites within Locations
##'   sc_transect_sigma = 0.02, # variability in Transects within Sites
##'   sc_sigma = 0.01, # random noise
##'   ma_site_sigma = 0.5, # variability in Sites within Locations
##'   ma_transect_sigma = 0.2, # variability in Transects within Sites
##'   ma_sigma = 0.1 # random noise
##' )
##' benthos_fixed_locs_obs <- sampling_design_fine_scale_fixed(benthos_fixed_locs_sf, config)
##' config <- list(
##'   Depths = 2,
##'   Depth_effect_multiplier = 2,
##'   Number_of_transects_per_site = 5,
##'   Number_of_frames_per_transect = 100,
##'   Points_per_frame = 5
##' )
##' benthos_fixed_locs_points <- sampling_design_fine_scale_points(benthos_fixed_locs_obs, config)
##' @export
sampling_design_fine_scale_points <- function(data_fixed_locs_obs, config) {
  set.seed(config$seed)
  ## put on fold scale
  data_fixed_locs_obs <-
    data_fixed_locs_obs |>
    tidyr::crossing(Depth = seq(3, 10, length = config$Depths)) |>
    tidyr::pivot_longer(cols = c(HCC, SC, MA),
      names_to = "Group",
      values_to = "Value") |>
    dplyr::group_by(Reef, Site, Transect, Year, Date) |>
    dplyr::mutate(Value = Value + rev(sort(config$Depth_effect_multiplier *
                                      scale(rnorm(config$Depths))))) |>
    dplyr::ungroup()
  ## Need to split the percentage cover into point and frames
  data_fixed_locs_obs <- data_fixed_locs_obs |>
    dplyr::group_by(Reef, Site, Transect, Year, Depth, Date) |>
    dplyr::mutate(
      Points = round(config$Number_of_frames_per_transect *
                       config$Points_per_frame *
                       (Value / sum(Value)), 0),
      Points = ifelse(Points < 0, 0, Points)
    ) |>
    tidyr::uncount(Points) |>
    dplyr::sample_n(dplyr::n(), replace = FALSE) |>
    dplyr::mutate(
      POINT_NO = rep_len(1:config$Points_per_frame, length = dplyr::n()),
      FRAME = rep(1:config$Number_of_frames_per_transect, each = config$Points_per_frame, length = dplyr::n())
    ) |>
    dplyr::ungroup()
  return(data_fixed_locs_obs)
}
