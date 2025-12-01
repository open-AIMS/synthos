#' Generate Benthic Cover Trajectory or Heatmap Plots
#'
#' Produces a set of plots showing temporal patterns in benthic cover across
#' survey depths and benthic classification groups. The function supports two
#' data formats: point-level benthic classifications (`"points"`) and
#' aggregated benthic cover values (`"cover"`).  
#' Depending on the `type` argument, the output is either a series of
#' time‐series trajectory plots or benthic‐cover heatmaps.
#'
#' @title Plot Benthic Cover Dynamics
#'
#' @param synthos_data A data frame containing benthic survey observations.
#'   Required variables depend on `data_type`:
#'   \itemize{
#'     \item `survey_depth` – depth category
#'     \item `project_name` – project identifier
#'     \item `site_name` – site name (e.g., "Reef10 Site 1")
#'     \item `survey_start_date` – survey timestamp
#'     \item `point_machine_classification` – benthic class (HCC, SC, MA, etc.)  
#'        *used when `data_type = "points"`*
#'     \item `cover` – numeric benthic cover values  
#'        *used when `data_type = "cover"`*
#'   }
#'
#' @param data_type A character string indicating the structure of input data:
#'   \itemize{
#'     \item `"points"` – raw point-based benthic classifications
#'     \item `"cover"` – pre-aggregated benthic cover values
#'   }
#'
#' @param type A character string specifying the type of plot to generate:
#'   \itemize{
#'     \item `"trajectories"` – time‐series plots (uses `plot_group()`)
#'     \item `"heatmaps"` – benthic cover heatmaps (uses `plot_map()`)
#'   }
#'
#' @return A named list of `ggplot` objects, where each plot corresponds to a
#'   unique combination of depth and benthic classification group. The resulting
#'   list can be visualised directly or iterated over using `purrr::walk()`.
#'
#' @details
#' Based on `data_type`, the function:
#'   \itemize{
#'     \item aggregates point classifications into benthic cover proportions; or
#'     \item averages numeric benthic cover values
#'   }
#'
#' Afterwards, it:
#'   \itemize{
#'     \item extracts year, reef, and site identifiers
#'     \item splits the dataset by survey depth and classification
#'     \item passes each subset to either `plot_group()` or `plot_map()`
#'     \item returns all plots as a list
#'   }
#'
#' @seealso [plot_group()], [plot_map()]
#' @author Julie
#' @export
plot_synthos <- function(synthos_data, type = "") {

  # ---- data prep: POINT-BASED DATA ----
  if (data_type == "points") {
    synthos_plot <- synthos_data |>
      dplyr::filter(!is.na(COUNT)) |>
      dplyr::group_by(
        survey_depth, project_name, reef, site,
        year, point_machine_classification
      ) |>
      dplyr::summarise(COUNT_site = sum(COUNT), TOTAL_site = sum(TOTAL) ,.groups = "drop_last") |>
      dplyr::mutate(
        COVER_site = (COUNT_site / TOTAL_site) * 100
      )
      # ) |>
      # dplyr::ungroup() |>
      # dplyr::mutate(
      #   year = lubridate::year(lubridate::ymd_hms(survey_start_date)),
      #   reef = stringr::str_extract(site_name, "^Reef\\d+"),
      #   site = stringr::str_extract(site_name, "Site \\d+$") |> stringr::str_remove("Site ")
      # )
  }

  # ---- data prep: AGGREGATED COVER ----
  if (data_type == "cover") {
    synthos_plot <- synthos_data |>
      dplyr::filter(!is.na(COVER)) |>
      dplyr::group_by(
        survey_depth, project_name, reef, site,
        year, point_machine_classification
      ) |>
      dplyr::summarise(COVER_site = mean(COVER), .groups = "drop") #|>
      # dplyr::mutate(
      #   year = lubridate::year(lubridate::ymd_hms(survey_start_date)),
      #   reef = stringr::str_extract(site_name, "^Reef\\d+"),
      #   site = stringr::str_extract(site_name, "Site \\d+$") |> stringr::str_remove("Site ")
      # )
  }

  # ---- common split variable ----
  split_group <- synthos_plot$point_machine_classification
  split_depth <- synthos_plot$survey_depth

  # ---- generate plots ----
  if (type == "trajectories") {
    plots_by_group <- synthos_plot |>
      split(list(split_depth, split_group)) |>
      purrr::map(plot_group)

  } else if (type == "heatmaps") {
    plots_by_group <- synthos_plot |>
      split(list(split_depth, split_group)) |>
      purrr::map(plot_map)
  }

  return(plots_by_group)
}
