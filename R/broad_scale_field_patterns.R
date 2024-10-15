baseline_hard_coral_cover <- function(spatial_grid, spde) {
  spatial_grid_pts_df <- spatial_grid |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    arrange(Longitude, Latitude)

  baseline_sample_hcc <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = clong + sin(clat) + # rnorm(1,0,1) +
        1.5 * clong + clat
    ) |>
    mutate(Y = scales::rescale(Y, to = c(-2, 0.8)))

  baseline_effects_hcc <- baseline_sample_hcc |>
    dplyr::select(Y) |>
    as.matrix()
  baseline_pts_sample_hcc <- inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    baseline_effects_hcc
  )
  baseline_pts_effects_hcc <- baseline_pts_sample_hcc |>
    cbind() |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    mutate(Year = as.numeric(Year))

  list(baseline_sample_hcc = baseline_sample_hcc,
    baseline_effects_hcc = baseline_effects_hcc,
    baseline_pts_sample_hcc = baseline_pts_sample_hcc,
    baseline_pts_effects_hcc = baseline_pts_effects_hcc)
}

synthetic_field_hcc <- function(all_effects_df, baseline_sample_hcc, spde) {
  spatial_grid_pts_df <- spatial_grid |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    arrange(Longitude, Latitude)
  ## Do all this on the link scale so that can use cumsum
  all_effects_hcc <- all_effects_df |>
    full_join(baseline_sample_hcc |>
                dplyr::select(Longitude, Latitude, BASE_HCC = Y)) |>
    group_by(Longitude, Latitude) |>
    mutate(HCC = BASE_HCC + Y_HCC) |>
    ungroup() |>
    dplyr::select(-BASE_HCC, -Y_HCC) |>
    pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = HCC
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  all_pts_sample_hcc <- inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects_hcc
  )
  all_pts_effects_hcc <- all_pts_sample_hcc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    mutate(
      Year = as.numeric(Year),
      Value = Value
    )
  list(all_effects_hcc = all_effects_hcc,
    all_pts_sample_hcc = all_pts_sample_hcc,
    all_pts_effects_hcc = all_pts_effects_hcc)
}


baseline_soft_coral_cover <- function(spatial_grid, spde) {
  spatial_grid_pts_df <- spatial_grid |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    arrange(Longitude, Latitude)

  baseline_sample_sc <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = clong + sin(clat) + # rnorm(1,0,1) +
        1.5 * clong + -1.5 * clat
    ) |>
    mutate(Y = scales::rescale(Y, to = c(-4, -2)))

  baseline_effects_sc <- baseline_sample_sc |>
    dplyr::select(Y) |>
    as.matrix()
  baseline_pts_sample_sc <- inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    baseline_effects_sc
  )
  baseline_pts_effects_sc <- baseline_pts_sample_sc |>
    cbind() |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    mutate(Year = as.numeric(Year))

  list(baseline_sample_sc = baseline_sample_sc,
    baseline_effects_sc = baseline_effects_sc,
    baseline_pts_sample_sc = baseline_pts_sample_sc,
    baseline_pts_effects_sc = baseline_pts_effects_sc)
}

synthetic_field_sc <- function(all_effects_df, baseline_sample_sc, spde) {
  spatial_grid_pts_df <- spatial_grid |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    arrange(Longitude, Latitude)
  ## Do all this on the link scale so that can use cumsum
  all_effects_sc <- all_effects_df |>
    full_join(baseline_sample_sc |>
                dplyr::select(Longitude, Latitude, BASE_SC = Y)) |>
    group_by(Longitude, Latitude) |>
    mutate(SC = BASE_SC + Y_SC) |>
    ungroup() |>
    dplyr::select(-BASE_SC, -Y_SC) |>
    pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = SC
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  all_pts_sample_sc <- inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects_sc
  )
  all_pts_effects_sc <- all_pts_sample_sc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    mutate(
      Year = as.numeric(Year),
      Value = Value
    )
  list(all_effects_sc = all_effects_sc,
    all_pts_sample_sc = all_pts_sample_sc,
    all_pts_effects_sc = all_pts_effects_sc)
}


synthetic_field_ma <- function(all_pts_effects_hcc, all_pts_effects_sc) {
  ## Do all this on the link scale so that can use cumsum
  all_pts_effects_ma <- all_pts_effects_hcc |>
    dplyr::rename(HCC=Value) |> 
    bind_cols(all_pts_effects_sc |>
                dplyr::select(SC=Value)) |>
    mutate(Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
      ## MA = Total_Avail*rbeta(n(), 2, 1),
      MA = Total_Avail,
      Value = qlogis(MA)) |>
    dplyr::select(-HCC, -SC, -Total_Avail, -MA)
  list(all_pts_effects_ma = all_pts_effects_ma)
}
