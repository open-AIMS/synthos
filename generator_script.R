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
surveys <- "random" # or "fixed"
data_type <- "points" #or "cover"

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

######################## Plots - cover at site level

plots <- plot_traj(synthos_data)
purrr::walk(plots, print)




# reef_data.synthetic_fixed_ready_site <- reef_data.synthetic_fixed_ready %>%
#  filter(!is.na(COUNT)) %>% filter(fGROUP == "HCC") %>% 
#  group_by(REEF_NAME, SITE_NO, fYEAR, fDEPTH, fGROUP) %>%
#  summarize(COUNT_sum = sum(COUNT),
#            TOTAL_sum = sum(TOTAL)) %>%
#  mutate(COVER_site = COUNT_sum / TOTAL_sum) %>%
#   mutate(fYEAR = as.numeric(as.character(fYEAR))) 

# if (length(unique(reef_data.synthetic_fixed_ready_site$fYEAR)) >10){
# p_heat <- ggplot(reef_data.synthetic_fixed_ready_site) +
#   geom_tile(aes(x = fYEAR, y = as.factor(SITE_NO),
#                 fill = COVER_site * 100)) +
#   scale_fill_viridis(
#     name = "Coral cover (%)", 
#     option = "plasma", 
#     begin = 0, 
#     end = ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE)), 
#     limits = c(0, ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE) * 100)), 
#     na.value = "grey90"
#   ) +
#   scale_x_continuous(
#     breaks = seq(min(reef_data.synthetic_fixed_ready_site$fYEAR, na.rm = TRUE),
#                  max(reef_data.synthetic_fixed_ready_site$fYEAR, na.rm = TRUE),
#                  by = 4)
#   ) +
#   facet_wrap(~REEF_NAME, ncol = 4) +
#   labs(x = "Year", y = "Site") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
#          axis.text.y = element_text(size=10),
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 10, margin = margin()),
#         legend.position = "bottom")
# }else{
# p_heat <- ggplot(reef_data.synthetic_fixed_ready_site) +
#   geom_tile(aes(x = fYEAR, y = as.factor(SITE_NO),
#                 fill = COVER_site * 100)) +
#   scale_fill_viridis(
#     name = "Coral cover (%)", 
#     option = "plasma", 
#     begin = 0, 
#     end = ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE)), 
#     limits = c(0, ceiling(max(reef_data.synthetic_fixed_ready_site$COVER_site, na.rm = TRUE) * 100)), 
#     na.value = "grey90"
#   ) +
#   facet_wrap(~REEF_NAME, ncol = 4) +
#   labs(x = "Year", y = "Site") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
#          axis.text.y = element_text(size=10),
#         strip.background = element_rect(fill = 'white'),
#         strip.text = element_text(size = 10, margin = margin()),
#         legend.position = "bottom")
# }

# ggsave(filename = paste0(title_of_run,"/report/extra/tile_data_",surveys,".png"),
#        plot = p_heat, width=13, height=12)  

