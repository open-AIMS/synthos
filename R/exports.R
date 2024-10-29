##' Prepare data for ReefCloud
##'
##' Prepare data for ReefCloud
##' @title Prepare data for ReefCloud
##' @param data_fixed_locs_points 
##' A sf object presenting the data of fixed locations and points
##' @return 
##' A data frame mimicking the export structure of reefCloud
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
##' reefcloud_synthetic_fixed_benthos <- prepare_for_reefcloud(benthos_fixed_locs_points)
##' @export
prepare_for_reefcloud <- function(data_fixed_locs_points) {
  reef_data_synthetic_fixed <-
    data_fixed_locs_points |>
    dplyr::mutate(
      project_id = 1,
      project_name = "synthetic_fixed",
      SITE_NO = stringr::str_replace(Site, "^S", "Site "),
      TRANSECT_NO = stringr::str_replace(Transect, "^T", "Transect "),
      site_name = factor(paste(Reef, SITE_NO)),
      site_id = as.numeric(site_name),
      site_latitude = Latitude,
      site_longitude = Longitude,
      site_depth = Depth,
      site_country = "synthetic Country",
      site_reef_name = factor(Reef),
      site_reef_type = NA,
      site_reef_zone = NA,
      site_code = NA,
      site_management = NA,
      survey_title = factor(paste(Reef, SITE_NO, TRANSECT_NO, format(Date, "%Y-%m-%d"))),
      survey_id = as.numeric(survey_title),
      survey_start_date = Date,
      survey_depth = Depth,
      survey_transect_number = as.numeric(stringr::str_replace(TRANSECT_NO, "Transect ", "")),
      ) 
  ## Photo-transect specific
  if ("POINT_NO" %in% names(reef_data_synthetic_fixed)) {
    reef_data_synthetic_fixed <-
      reef_data_synthetic_fixed |>
      mutate(
        image_name = factor(paste(survey_title, FRAME)),
        image_id = as.numeric(image_name),
        image_quality = 100,
        point_no = POINT_NO,
        point_id = as.numeric(factor(paste(image_name, POINT_NO))),
        point_machine_classification = Group
      )
  }
  ## Quadrat (%cover) specific
  if ("QUAD_NO" %in% names(reef_data_synthetic_fixed)) {
    reef_data_synthetic_fixed <-
      reef_data_synthetic_fixed |> 
      mutate(
        cover = Value,
        quad_no = as.numeric(factor(QUAD_NO))
      )
  }
  reef_data_synthetic_fixed <-
    reef_data_synthetic_fixed |> 
        dplyr::select(
          project_id,
          project_name,
          site_id,
          site_name,
          site_latitude,
          site_longitude,
          site_depth,
          site_country,
          site_reef_name,
          site_reef_type,
          site_reef_zone,
          site_code,
          site_management,
          survey_id,
          survey_title,
          survey_start_date,
          survey_depth,
          survey_transect_number,
          any_of(c(
            "image_id",
            "image_name",
            "image_quality",
            "point_id",
            "point_no",
            "point_machine_classification"
          )),
          any_of(c(
            "quad_no",
            "cover"
          ))
        )
  return(reef_data_synthetic_fixed)
}
