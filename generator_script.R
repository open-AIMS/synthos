## Generating synthethic data
rm(list = ls())


library(sf)
library(stars)
library(gstat)
library(INLA)
library(ggplot2)
detach("package:synthos", unload = TRUE, character.only = TRUE)
remotes::install_github("open-AIMS/synthos@julie", force = TRUE, dependencies = FALSE)
library(synthos)

##### Generate settings
surveys <- "random"
data_type <- "points"

synthos::generateSettings(nreefs = 25, nsites = 3, nyears = 15)

# 1. Generation of the spatial and temporal domains, disturbances and baselines
spatial_domain <- st_geometry(
  st_multipoint(
    x = rbind(
      c(0, -11),
      c(3,-11),
      c(6,-14),
      c(1,-15),
      c(2,-12),
      c(0,-11)
    )
  )
) |>
  st_set_crs(config_sp$crs) |>
  st_cast("POLYGON")

## ---- SpatialPoints
set.seed(config_sp$seed)
spatial_grid <- spatial_domain |>
  st_set_crs(NA) |>
  st_sample(size = 10000, type = "regular") |>
  st_set_crs(config_sp$crs)
sf_use_s2(FALSE)

benthos_reefs_pts <- synthos::create_synthetic_reef_landscape(spatial_grid, config_sp)
## ----end


##-----------------------------#
## 1. Select sampling design
##-----------------------------#

if (surveys == "fixed") {
  locs_sf <- synthos::sampling_design_large_scale_fixed(
    benthos_reefs_pts, config_lrge
  )
  obs <- synthos::sampling_design_fine_scale_fixed(
    locs_sf, config_fine
  )
  
} else if (surveys == "random") {
  locs_sf <- synthos::sampling_design_large_scale_random(
    benthos_reefs_pts, config_lrge
  )
  obs <- synthos::sampling_design_fine_scale_random(
    locs_sf, config_fine
  )
} else {
  stop("surveys must be 'fixed' or 'random'.")
}

##-----------------------------#
## 2. Generate export table
##-----------------------------#

if (data_type == "points") {
  pts <- synthos::sampling_design_fine_scale_points(obs, config_pt)
  synthos_data <- synthos::prepare_table(pts)
  
} else if (data_type == "cover") {
  cov <- synthos::sampling_design_fine_scale_cover(obs, config_pt)
  synthos_data <- synthos::prepare_table(cov)
  
} else {
  stop("data_type must be 'points' or 'cover'.")
}

######################## Plots - cover at transect level

if (data_type == "points") {
synthos_plot <- synthos_data |>
  dplyr::group_by(survey_depth, project_name, site_name, survey_transect_number, survey_start_date, point_machine_classification) |>
  dplyr::summarise(COUNT = dplyr::n()) |>
  dplyr::ungroup(point_machine_classification) |>
  dplyr::mutate(TOTAL=sum(COUNT)) |>
  dplyr::ungroup() |>
  dplyr::mutate(COVER = COUNT / TOTAL) |>
  dplyr::mutate(year = lubridate::year(lubridate::ymd_hms(survey_start_date))) |>
  dplyr::mutate(reef = stringr::str_extract(site_name, "^Reef\\d+")) |>
  dplyr::mutate(site = stringr::str_extract(site_name, "Site \\d+$") |> stringr::str_remove("Site "))

ggplot2::ggplot(synthos_data |> dplyr::filter(survey_depth == "10") |> dplyr::filter(point_machine_classification == "HCC")) + 
  ggplot2::geom_line(ggplot2::aes(x = year, y = COVER*100, group = interaction(as.factor(survey_transect_number), as.factor(reef), as.factor(site)),
  col = as.factor(site_name)), 
   show.legend = FALSE) + 
  ggplot2::facet_wrap(~reef, ncol=4) + ggplot2::theme_bw() +
  ggplot2::labs(x = "Year", y = "Coral cover") +
  ggplot2::ylab("Coral cover") + ggplot2::xlab("Year") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = ggplot2::element_text(size=10),
        axis.title.y = ggplot2::element_text(size=11),
        axis.title.x= ggplot2::element_text(size=11),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_rect(fill = 'white'),
        strip.text = ggplot2::element_text(size = 10, margin = ggplot2::margin())) + 
  ggplot2::ggtitle("Fixed design")

}


## Data viz 

p_vis_data_fixed <- ggplot(reef_data.synthetic_fixed_ready %>% filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC")) + 
  geom_line(aes(x = fYEAR, y = COVER*100, group = interaction(as.factor(TRANSECT_NO), as.factor(SITE_NO), REEF_NAME),
  col = as.factor(SITE_NO)), 
   show.legend = FALSE) + 
  facet_wrap(~REEF_NAME, ncol=4) + theme_bw() +
  labs(x = "Year", y = "Coral cover") +
  ylab("Coral cover") + xlab("Year")+theme_bw()+
  theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1),legend.position = "right",
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(size=11),
        axis.title.x= element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin())) + 
  ggtitle("Fixed design")

ggsave(filename = paste0(title_of_run,"/report/extra/trend_data_",surveys,".png"),
       plot = p_vis_data_fixed, width=13, height=12)  

reef_data.synthetic_fixed_ready_site <- reef_data.synthetic_fixed_ready %>%
 filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC") %>% 
 group_by(REEF_NAME, SITE_NO, fYEAR, fDEPTH, fGROUP) %>%
 summarize(COUNT_sum = sum(COUNT),
           TOTAL_sum = sum(TOTAL)) %>%
 mutate(COVER_site = COUNT_sum / TOTAL_sum) %>%
  mutate(fYEAR = as.numeric(as.character(fYEAR))) 

if (length(unique(reef_data.synthetic_fixed_ready_site$fYEAR)) >10){
p_heat <- ggplot(reef_data.synthetic_fixed_ready_site) +
  geom_tile(aes(x = fYEAR, y = as.factor(SITE_NO),
                fill = COVER_site * 100)) +
  scale_fill_viridis(
    name = "Coral cover (%)", 
    option = "plasma", 
    begin = 0, 
    end = ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE)), 
    limits = c(0, ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE) * 100)), 
    na.value = "grey90"
  ) +
  scale_x_continuous(
    breaks = seq(min(reef_data.synthetic_fixed_ready_site$fYEAR, na.rm = TRUE),
                 max(reef_data.synthetic_fixed_ready_site$fYEAR, na.rm = TRUE),
                 by = 4)
  ) +
  facet_wrap(~REEF_NAME, ncol = 4) +
  labs(x = "Year", y = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
         axis.text.y = element_text(size=10),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin()),
        legend.position = "bottom")
}else{
p_heat <- ggplot(reef_data.synthetic_fixed_ready_site) +
  geom_tile(aes(x = fYEAR, y = as.factor(SITE_NO),
                fill = COVER_site * 100)) +
  scale_fill_viridis(
    name = "Coral cover (%)", 
    option = "plasma", 
    begin = 0, 
    end = ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE)), 
    limits = c(0, ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE) * 100)), 
    na.value = "grey90"
  ) +
  facet_wrap(~REEF_NAME, ncol = 4) +
  labs(x = "Year", y = "Site") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
         axis.text.y = element_text(size=10),
        strip.background = element_rect(fill = 'white'),
        strip.text = element_text(size = 10, margin = margin()),
        legend.position = "bottom")
}

ggsave(filename = paste0(title_of_run,"/report/extra/tile_data_",surveys,".png"),
       plot = p_heat, width=13, height=12)  

