#' Generate Trajectory Plots for Benthic Cover
#'
#' Creates a list of trajectory plots showing temporal changes in benthic cover
#' across survey depths and benthic classifications. The function handles two
#' data types: point-based benthic classifications ("points") and aggregated
#' percentage cover ("cover"). It returns a nested list of `ggplot` objects,
#' where each plot corresponds to a unique combination of depth and
#' point-machine-classification.
#'
#' @title Plot Benthic Cover Trajectories
#'
#' @param synthos_data A data frame containing benthic survey observations. Must
#'   include (depending on `data_type`):
#'   \itemize{
#'     \item `survey_depth` – depth category
#'     \item `project_name` – project identifier
#'     \item `site_name` – site identifier (e.g., "Reef10 Site 1")
#'     \item `survey_start_date` – timestamp of survey
#'     \item `point_machine_classification` – benthic class (HCC, SC, MA, etc.) for `data_type = "points"`
#'     \item `cover` – percentage cover values for `data_type = "cover"`
#'   }
#'
#' @param data_type A string specifying the data source structure:
#'   \itemize{
#'     \item `"points"` – uses raw point-level machine classifications
#'     \item `"cover"` – uses pre-aggregated percentage cover values
#'   }
#'
#' @return A named list of `ggplot` objects, each representing a benthic cover
#'   trajectory over time for a unique combination of depth and classification
#'   group. Plots can be viewed individually or iterated over using
#'   `purrr::walk()`.
#'
#' @details
#' Depending on `data_type`, this function:
#'   \itemize{
#'     \item aggregates point-based classifications into cover proportions, or
#'     \item averages numeric cover values
#'   }
#'
#' It then:
#'   \itemize{
#'     \item extracts temporal and site identifiers (year, reef, site)
#'     \item splits the dataset by depth and classification
#'     \item applies `plot_group()` to each data subset
#'     \item returns all plots as a list
#'   }
#'
#' @seealso [plot_group()]
#' @author Julie
##' @export
plot_traj <- function(synthos_data){

  if (data_type == "points") {

    synthos_plot <- synthos_data |>
      dplyr::group_by(survey_depth, project_name, site_name, survey_start_date, point_machine_classification) |>
      dplyr::summarise(COUNT = dplyr::n()) |>
      dplyr::ungroup(point_machine_classification) |>
      dplyr::mutate(TOTAL = sum(COUNT)) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        COVER = (COUNT / TOTAL) * 100,
        year  = lubridate::year(lubridate::ymd_hms(survey_start_date)),
        reef  = stringr::str_extract(site_name, "^Reef\\d+"),
        site  = stringr::str_extract(site_name, "Site \\d+$") |> stringr::str_remove("Site ")
      )


  } else if (data_type == "cover") {

    synthos_plot <- synthos_data |>
      dplyr::group_by(survey_depth, project_name, site_name, survey_start_date, point_machine_classification) |>
      dplyr::summarise(COVER = mean(cover)) |>
      dplyr::mutate(
        year = lubridate::year(lubridate::ymd_hms(survey_start_date)),
        reef = stringr::str_extract(site_name, "^Reef\\d+"),
        site = stringr::str_extract(site_name, "Site \\d+$") |> stringr::str_remove("Site ")
      )
  }

  # ---- final split and plotting ----
    split_depth  <- synthos_plot$survey_depth
    split_group  <- synthos_plot$point_machine_classification

  plots_by_group <- synthos_plot %>% 
    split(list(split_depth, split_group)) %>%
    purrr::map(plot_group)

  return(plots_by_group)
}

