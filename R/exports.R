
prepare_for_reefcloud <- function(data_fixed_locs_obs) {
  reef_data_synthetic_fixed <-
    data_fixed_locs_obs |>
    mutate(
      project_id = 1,
      project_name = "synthetic_fixed",
      SITE_NO = str_replace(Site, "^S", "Site "),
      TRANSECT_NO = str_replace(Transect, "^T", "Transect "),
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
      survey_transect_number = as.numeric(str_replace(TRANSECT_NO, "Transect ", "")),
      image_name = factor(paste(survey_title, FRAME)),
      image_id = as.numeric(image_name),
      image_quality = 100,
      point_no = POINT_NO,
      point_id = as.numeric(factor(paste(image_name, POINT_NO))),
      point_machine_classification = Group
    ) |>
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
      image_id,
      image_name,
      image_quality,
      point_id,
      point_no,
      point_machine_classification
    )
  return(reef_data_synthetic_fixed)
}
