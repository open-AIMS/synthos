##' Large scale fixed sampling design
##'
##' Given a number of locations and a number of sites per location,
##' select the locations from the full field (domain) and then
##' randomly select the sites within those locations.
##'
##' @title Large scale fixed sampling design
##' @param data_reefs_pts_sf
##' A sf containing the full field with:
##' - Year: a numeric representation of years in the temporal domain
##' - Reef: a unique reef id
##' - HCC: hard coral cover (on logit scale)
##' - SC: soft coral cover (on logit scale)
##' - MA: macroalgae cover (on logit scale)
##' - geometry: sf geometry
##' @param config
##' A list containing:
##' - n_locs: number of reefs (locations)
##' - n_sites: number of sites per location
##' - seed: random seed
##' @return
##' A sf object denoting the large scale sampling design.  Cover values are on a logit scale.
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
##'   hcc_growth = 0.3,
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
##' @export
sampling_design_large_scale_fixed <- function(data_reefs_pts_sf, config) {
  testthat::expect(
    inherits(data_reefs_pts_sf, c("sf")),
    "data_reefs_pts_sf must be an sf object"
  )
  testthat::expect_in(
    sort(c(
      "seed",
      "n_locs",
      "n_sites"
    )),
    sort(names(config))
  )
  set.seed(config$seed)
  ## Start by randomly selecting nLocs Reefs
  Reefs_fixed <- data_reefs_pts_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(Reef) |>
    dplyr::distinct() |>
    dplyr::sample_n(size = config$n_locs) |>
    dplyr::pull(Reef)
  ## Then filter to these Reefs before selecting a single location within
  ## each of the Reefs
  data_fixed_locs_sf <- data_reefs_pts_sf |>
    dplyr::filter(Reef %in% Reefs_fixed) |>
    dplyr::select(Reef, geometry) |>
    dplyr::distinct(.keep_all = TRUE) |>
    dplyr::group_by(Reef) |>
    dplyr::sample_n(config$n_sites) |>
    dplyr::mutate(Site = paste0("S", 1:dplyr::n())) |>
    dplyr::ungroup() |>
    sf::st_join(data_reefs_pts_sf |>
      dplyr::select(-Reef))
  return(data_fixed_locs_sf)
}

##' Fine scale fixed sampling design
##'
##' Given a large scale sampling design and a number of parametes,
##' create a spatial hierarchy of transects per site, photo-frames
##' per transect and points per frame.
##' @title  Fine scale fixed sampling design
##' @param data_fixed_locs_sf
##' A sf containing the large scale sampling design with:
##' - Reef: a unique reef id
##' - Site: a unique site id
##' - geometry: sf geometry
##' @param config
##' A list containing:
##' - Number_of_transects_per_site: number of transects per site
##' - Depths: number of depths
##' - Number_of_frames_per_transect: number of frames per transect
##' - Points_per_frame: number of points per frame
##' - hcc_site_sigma: variability (on logit scale) in Sites within Locations
##' - hcc_transect_sigma: variability (on logit scale) in Transects within Sites
##' - hcc_sigma: random noise (on a logit scale)
##' - sc_site_sigma: variability (on a logit scale) in Sites within Locations
##' - sc_transect_sigma: variability (on a logit scale) in Transects within Sites
##' - sc_sigma: random noise (on a logit scale)
##' - ma_site_sigma: variability (on a logit scale) in Sites within Locations
##' - ma_transect_sigma: variability (on a logit scale in Transects within Sites
##' - ma_sigma: random noise (on a logit scale
##' @return
##' A sf object denoting the fine scale sampling design.  Cover values are on natural percentage cover scale
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
##' @export
sampling_design_fine_scale_fixed <- function(data_fixed_locs_sf, config) {
  set.seed(config$seed)
  data_fixed_locs_obs <- data_fixed_locs_sf |>
    dplyr::bind_cols(data_fixed_locs_sf |>
                sf::st_coordinates() |>
                as.data.frame() |>
                dplyr::rename(Longitude = X, Latitude = Y)) |>
    sf::st_drop_geometry() |>
    as.data.frame() |>
    dplyr::group_by(Longitude, Latitude, Reef) |>
    tidyr::crossing(
      Transect = paste0("T",1:config$Number_of_transects_per_site)) |>
    dplyr::group_by(Site, .add = TRUE) |>
    dplyr::mutate(
      SiteEffects_HCC = rnorm(1, 0, config$hcc_site_sigma),
      SiteEffects_SC = rnorm(1, 0, config$sc_site_sigma),
      SiteEffects_MA = rnorm(1, 0, config$ma_site_sigma)
    ) |>
    dplyr::group_by(Transect, .add = TRUE) |>
    dplyr::mutate(
      TransectEffects_HCC = rnorm(1, 0, config$hcc_transect_sigma),
      TransectEffects_SC = rnorm(1, 0, config$sc_transect_sigma),
      TransectEffects_MA = rnorm(1, 0, config$ma_transect_sigma)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      HCC1 = HCC + SiteEffects_HCC +
        TransectEffects_HCC +
        rnorm(dplyr::n(), 0, config$hcc_sigma),
      HCC2 = 100*plogis(HCC1),
      SC1 = SC + SiteEffects_SC + TransectEffects_SC +
        rnorm(dplyr::n(), 0, config$sc_sigma),
      SC2 = 100*plogis(SC1),
      MA1 = MA + SiteEffects_MA + TransectEffects_MA
      + rnorm(dplyr::n(), 0, config$ma_sigma),
      MA2 = 100*plogis(MA1)
    ) |>
    dplyr::arrange(Reef, Site, Transect, Year) |>
    dplyr::select(Reef, Longitude, Latitude, Site,
      Transect, Year, HCC = HCC2, SC = SC2, MA = MA2) |>
    ## dplyr::mutate(Year = 2021 - max(config$years) + Year,
    dplyr::mutate(Date = as.POSIXct(paste0(Year, "-01-01 14:00:00")))
  return(data_fixed_locs_obs)
}

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

##' Fine scale % cover estimates
##'
##' Disperse benthos cover across quadrats
##' @title Fine scale percent cover
##' @param data_fixed_locs_obs 
##' A sf object representing ...
##' @param config
##' A list containing:
##' - Depths: number of depths
##' - Depth_effect_multiplier: a multiplier for the depth effect
##' - Number_of_transects_per_site: number of transects per site
##' - Number_of_quadrats_per_transect: number of quadrats per transect
##' - Quad_sigma: standard deviation of the quadrat effect
##' @return 
##' A sf object representing the fine scale sampling design with quadrats and percentage cover
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
##'   Quad_sigma = 1
##' )
##' benthos_fixed_locs_points <- sampling_design_fine_scale_points(benthos_fixed_locs_obs, config)
##' @export
sampling_design_fine_scale_cover <- function(data_fixed_locs_obs, config) {
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
  ## Need to split the percentage cover into quadrats
  data_fixed_locs_obs <- data_fixed_locs_obs |>
    dplyr::group_by(Reef, Site, Transect, Year, Depth, Date) |>
    tidyr::crossing(
      Quad = paste0("Q", 1:config$Number_of_quadrats_per_transect)
    ) |>
    mutate(Value = 100* plogis(qlogis(Value / 100) + rnorm(n(), 0, config$Quad_sigma))) |> 
    dplyr::ungroup()
  return(data_fixed_locs_obs)
}

