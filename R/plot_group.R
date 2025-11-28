#' Plot Coral Cover Time Series by Group
#'
#' Generates a faceted time‐series plot of coral cover (%) for a given subset
#' of the dataset. The plot shows temporal dynamics across sites, grouped by
#' machine classification (e.g., HCC, SC, MA) and survey depth. A dynamic
#' colour palette is used so that site colours adjust automatically to the
#' number of unique sites present in the data slice.
#'
#' @title Plot Time Series of Coral Cover by Group and Depth
#'
#' @param df A data frame containing at least the following variables:
#'   \itemize{
#'     \item `year` – survey year (numeric)
#'     \item `COVER` – coral cover proportion (0–1)
#'     \item `reef` – reef identifier
#'     \item `site` – site identifier
#'     \item `survey_depth` – depth category (character or numeric)
#'     \item `point_machine_classification` – benthic classification label
#'   }
#'
#' @return A `ggplot` object showing coral cover over time, faceted by reef
#'   and coloured by site.
#'
#' @details
#' This function:
#'   \itemize{
#'     \item creates a dynamic site-level colour palette using `scales::hue_pal()`
#'     \item plots coral cover time series at the site level
#'     \item facets results across reefs for visual comparison
#'     \item automatically labels panels with depth and classification group
#'   }
#'
#' @author Julie
#' @export
#' 
plot_group <- function(df) {

  sites <- sort(unique(df$site))
  n_sites <- length(sites)

  ggplot(df) + 
    geom_line(aes(
      x = year, 
      y = COVER * 100, 
      group = site, 
      color = site
    )) +
    facet_wrap(~ reef, ncol = 5) +
    scale_color_manual(values = scales::hue_pal(l = 30)(n_sites)) +
    theme_bw() +
    labs(
      x = "Year", 
      y = "Coral cover (%)", 
      color = "Site",
      title = paste0("Group: ", unique(df$point_machine_classification), " and ", "Depth: ", unique(df$survey_depth))
    ) +
    theme(
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      axis.title.x = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    )
}
