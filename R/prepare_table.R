#' Prepare Synthetic ReefCloud Data Table
#'
#' Converts point-level benthic sampling observations into a structured
#' data frame mimicking the ReefCloud pipeline export format
#' (https://reefcloud.ai/). Assigns survey, site, transect, frame, and point
#' identifiers, along with benthic cover values and relevant metadata.
#'
#' @title Prepare ReefCloud-Compatible Data Table
#'
#' @param data_fixed_locs_points An `sf` object containing fine-scale point-level
#' sampling data produced by `sampling_design_fine_scale_points()`. Must include:
#' \itemize{
#'   \item Reef, Site, Transect, Year, Date
#'   \item Depth, Frame, POINT_NO (point number)
#'   \item Benthic cover values and Group (HCC, SC, MA)
#'   \item Coordinates (geometry)
#' }
#'
#' @return A `data.frame` formatted for ReefCloud export, containing:
#' \itemize{
#'   \item Project, site, survey, transect identifiers
#'   \item Geographic coordinates and depth
#'   \item Image/frame identifiers (if photo-transect data)
#'   \item Point identifiers and benthic cover (Value or Group)
#'   \item Quadrat information if available
#' }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Generates consistent site, survey, and transect IDs
#'   \item Adds frame and point identifiers for photo-point data
#'   \item Adds cover and quadrat identifiers for quadrat data
#'   \item Selects only columns relevant to ReefCloud export
#' }
#'
#' @author Murray
#' @export
prepare_table <- function(data_fixed_locs_points) {
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
      dplyr::mutate(
        image_name = factor(paste(survey_title, FRAME)),
        image_id = as.numeric(image_name),
        image_quality = 100,
        point_no = POINT_NO,
        point_id = as.numeric(factor(paste(image_name, POINT_NO))),
        point_machine_classification = Group
      )
  }
  ## Quadrat (%cover) specific
  if ("Quad" %in% names(reef_data_synthetic_fixed)) {
    reef_data_synthetic_fixed <-
      reef_data_synthetic_fixed |> 
      dplyr::mutate(
        cover = Value,
        quad_no = as.numeric(factor(Quad)),
        point_machine_classification = Group
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
