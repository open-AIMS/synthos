sampling_design_large_scale <- function(data_reefs_pts_sf, config) {
  set.seed(config$seed)
  ## Start by randomly selecting nLocs Reefs
  Reefs_fixed <- data_reefs_pts_sf |>
    st_drop_geometry() |>
    dplyr::select(Reef) |>
    distinct() |>
    sample_n(size = config$n_locs) |>
    pull(Reef)
  ## Then filter to these Reefs before selecting a single location within
  ## each of the Reefs
  data_fixed_locs_sf <- data_reefs_pts_sf |>
    filter(Reef %in% Reefs_fixed) |>
    dplyr::select(Reef, geometry) |>
    distinct(.keep_all = TRUE) |> 
    group_by(Reef) |>
    sample_n(config$n_sites) |>
    mutate(Site = paste0("S", 1:n())) |>
    ungroup() |> 
    st_join(data_reefs_pts_sf |> 
              dplyr::select(-Reef))

  return(data_fixed_locs_sf)
}

sampling_design_fine_scale <- function(data_fixed_locs_sf, config) {
  set.seed(config$seed)
  data_fixed_locs_obs <- data_fixed_locs_sf |>
    bind_cols(data_fixed_locs_sf |>
                st_coordinates() |>
                as.data.frame() |>
                dplyr::rename(Longitude = X, Latitude = Y)) |>
    st_drop_geometry() |>
    as.data.frame() |>
    group_by(Longitude, Latitude, Reef) |>
    crossing(
      Transect = paste0("T",1:config$Number_of_transects_per_site)) |>
    group_by(Site, .add = TRUE) |>
    mutate(
      SiteEffects_HCC = rnorm(1, 0, config$hcc_site_sigma),
      SiteEffects_SC = rnorm(1, 0, config$sc_site_sigma),
      SiteEffects_MA = rnorm(1, 0, config$ma_site_sigma)
    ) |>
    group_by(Transect, .add = TRUE) |>
    mutate(
      TransectEffects_HCC = rnorm(1, 0, config$hcc_transect_sigma),
      TransectEffects_SC = rnorm(1, 0, config$sc_transect_sigma),
      TransectEffects_MA = rnorm(1, 0, config$ma_transect_sigma)
    ) |>
    ungroup() |>
    mutate(
      HCC1 = HCC + SiteEffects_HCC +
        TransectEffects_HCC +
        rnorm(n(), 0, config$hcc_sigma),
      HCC2 = 100*plogis(HCC1),
      SC1 = SC + SiteEffects_SC + TransectEffects_SC +
        rnorm(n(), 0, config$sc_sigma),
      SC2 = 100*plogis(SC1),
      MA1 = MA + SiteEffects_MA + TransectEffects_MA
      + rnorm(n(), 0, config$ma_sigma),
      MA2 = 100*plogis(MA1)
    ) |>
    arrange(Reef, Site, Transect, Year) |>
    dplyr::select(Reef, Longitude, Latitude, Site,
      Transect, Year, HCC = HCC2, SC = SC2, MA = MA2) |>
    mutate(Year = 2021 - max(config$years) + Year,
      Date = as.POSIXct(paste0(Year, "-01-01 14:00:00")))
  return(data_fixed_locs_obs)
}

sampling_design_fine_scale_points <- function(data_fixed_locs_obs, config) {
  set.seed(config$seed)
  ## put on fold scale
  data_fixed_locs_obs <-
    data_fixed_locs_obs |>
    tidyr::crossing(Depth = seq(3, 10, length = config$Depths)) |>
    pivot_longer(cols = c(HCC, SC, MA),
      names_to = "Group",
      values_to = "Value") |>
    group_by(Reef, Site, Transect, Year, Date) |>
    mutate(Value = Value + rev(sort(config$Depth_effect_multiplier *
                                      scale(rnorm(config$Depths))))) |>
    ungroup()
  ## Need to split the percentage cover into point and frames
  data_fixed_locs_obs <- data_fixed_locs_obs |>
    group_by(Reef, Site, Transect, Year, Depth, Date) |>
    mutate(
      Points = round(config$Number_of_frames_per_transect *
                       config$Points_per_frame *
                       (Value / sum(Value)), 0),
      Points = ifelse(Points < 0, 0, Points)
    ) |>
    tidyr::uncount(Points) |>
    sample_n(n(), replace = FALSE) |>
    mutate(
      POINT_NO = rep_len(1:config$Points_per_frame, length = n()),
      FRAME = rep(1:config$Number_of_frames_per_transect, each = config$Points_per_frame, length = n())
    ) |>
    ungroup()
  return(data_fixed_locs_obs)
}

