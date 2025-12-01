#' Plot Benthic COVER_site Heatmap by Reef and Site
#'
#' Creates a faceted heatmap showing benthic COVER_site (%) through time for a
#' selected subset of the data. Each tile represents the benthic COVER_site at a
#' given site and year. The plot is faceted by reef, with an automatic colour
#' scale and optional custom year breaks for long time series.
#'
#' @title Benthic COVER_site Heatmap
#'
#' @param df A data frame containing:
#'   \itemize{
#'     \item `year` – survey year (numeric)
#'     \item `COVER_site` – benthic COVER_site percentage (0–100)
#'     \item `reef` – reef identifier
#'     \item `site` – site identifier
#'   }
#'
#' @return A `ggplot` object visualising benthic COVER_site as a heatmap, faceted by
#'   reef.
#'
#' @details
#' The function:
#'   \itemize{
#'     \item uses a viridis colour scale for benthic COVER_site
#'     \item applies dynamic year breaks when the time series exceeds 10 years
#'     \item facets results by reef for comparison across locations
#'   }
#'
#' @author Julie
#' @export
plot_map <- function(df) {


if (length(unique(df$year)) >10){
  ggplot(df) +
  geom_tile(aes(x = year, y = as.factor(site),
                fill = COVER_site)) +
  viridis::scale_fill_viridis(
    name = "Cover (%)", 
    option = "plasma", 
    begin = 0, 
    end = ceiling(max(df$COVER_site/100, na.rm = TRUE)), 
    limits = c(0, ceiling(max(df$COVER_site, na.rm = TRUE))), 
    na.value = "grey90"
  ) +
  scale_x_continuous(
    breaks = seq(min(df$year, na.rm = TRUE),
                 max(df$year, na.rm = TRUE),
                 by = 4)
  ) +
  facet_wrap(~reef, ncol = 4) +
  labs(x = "Year", y = "",
      title = paste0("Group: ", unique(df$point_machine_classification), " and ", "Depth: ", unique(df$survey_depth))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
         axis.text.y = element_text(size=10),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin()),
        legend.position = "bottom")
}else{
  ggplot(df) +
  geom_tile(aes(x = year, y = as.factor(site),
                fill = COVER_site)) +
  viridis::scale_fill_viridis(
    name = "Cover (%)", 
    option = "plasma", 
    begin = 0, 
    end = ceiling(max(df$COVER_site/100, na.rm = TRUE)), 
    limits = c(0, ceiling(max(df$COVER_site, na.rm = TRUE))), 
    na.value = "grey90"
  ) +
  facet_wrap(~reef, ncol = 4) +
  labs(x = "Year", y = "",
      title = paste0("Group: ", unique(df$point_machine_classification), " and ", "Depth: ", unique(df$survey_depth))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
         axis.text.y = element_text(size=10),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin()),
        legend.position = "bottom")
}

}
