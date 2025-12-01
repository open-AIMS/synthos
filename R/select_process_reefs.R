#' Select Reefs and Prepare Synthetic Benthic Data
#'
#' Processes a synthetic reef dataset to select a subset of reefs for survey
#' and prepares benthic cover data either at the point level or aggregated cover
#' level. Reef and site identifiers are extracted, and site names are
#' standardized.
#'
#' @title Select and Process Synthetic Reefs
#'
#' @param reef_data_synthetic A data frame containing synthetic reef survey data.
#'   Must include columns such as `site_name`, `survey_depth`, `survey_start_date`,
#'   `survey_transect_number`, and `point_machine_classification`. If `data_type = "cover"`,
#'   it should also include `cover`.
#'
#' @param data_type A string specifying the type of data to process:
#'   \itemize{
#'     \item `"points"` – process raw point-based benthic classifications
#'     \item `"cover"` – process aggregated percentage cover values
#'   }
#'
#' @return A processed data frame containing:
#'   \itemize{
#'     \item `reef` – reef identifier
#'     \item `site` – standardized site name ("Site X")
#'     \item `year` – survey year extracted from `survey_start_date`
#'     \item `COUNT` – selected point counts (for points data)
#'     \item `COVER` – benthic cover (%) for selected reefs
#'   }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Randomly selects reefs based on `config_lrge$n_locs`
#'   \item Aggregates point-based counts or cover values depending on `data_type`
#'   \item Extracts year, reef, and site identifiers from `site_name`
#'   \item Standardizes site names from "S1" to "Site 1"
#'   \item Assigns NA to reefs not selected for survey
#' }
#'
#' @author Julie
#' @export
select_process_reefs <- function(reef_data_synthetic) {

  # ---- select reefs to survey ----
  set.seed(config_lrge$seed)

  reefs_selected <- reef_data_synthetic |>
    dplyr::mutate(reef = stringr::str_extract(site_name, "^Reef\\d+")) |>
    dplyr::select(reef) |>
    dplyr::distinct() |>
    dplyr::sample_n(size = config_lrge$n_locs) |>
    dplyr::pull(reef)

  # ---- data prep: POINT-BASED DATA ----
  if (data_type == "points") {
    reef_data_synthetic <- reef_data_synthetic |>
      dplyr::group_by(
        survey_depth, project_name, site_name, survey_transect_number,
        survey_start_date, point_machine_classification
      ) |>
      dplyr::summarise(
        COUNT_TRUE = dplyr::n(),
        .groups = "drop_last"
      ) |>
      dplyr::mutate(
        TOTAL = sum(COUNT_TRUE, na.rm = TRUE),
        COVER_TRUE = (COUNT_TRUE / TOTAL) * 100
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        year = lubridate::year(lubridate::ymd_hms(survey_start_date)),
        reef = stringr::str_extract(site_name, "^Reef\\d+"),
        site = stringr::str_extract(site_name, "S\\d+$") |>
               stringr::str_remove("S") |>
               (\(x) paste("Site", x))()
      ) |>
      dplyr::mutate(
        COUNT = dplyr::case_when(
          reef %in% reefs_selected ~ COUNT_TRUE,
          TRUE ~ NA_real_
        ),
        COVER = COUNT / TOTAL
      )
  }

  # ---- data prep: AGGREGATED COVER ----
  if (data_type == "cover") {
    reef_data_synthetic <- reef_data_synthetic |>
      dplyr::group_by(
        survey_depth, project_name, site_name, survey_transect_number,
        survey_start_date, point_machine_classification
      ) |>
      dplyr::summarise(COVER_TRUE = mean(cover, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(
        year = lubridate::year(lubridate::ymd_hms(survey_start_date)),
        reef = stringr::str_extract(site_name, "^Reef\\d+"),
        site = stringr::str_extract(site_name, "S\\d+$") |>
               stringr::str_remove("S") |>
               (\(x) paste("Site", x))()
      ) |>
      dplyr::mutate(
        COVER = dplyr::case_when(
          reef %in% reefs_selected ~ COVER_TRUE,
          TRUE ~ NA_real_
        )
      )
  }

  return(reef_data_synthetic)
}
