rasterize_reefs <- function(reefs_sf) {
  data_reefs_sf <- reefs_sf |>
    st_as_stars(dx = 0.01) |>  # rasterize
    st_as_sf(as_points = TRUE) |>
    filter(values == 1L)

  data_reefs_df <- data_reefs_sf |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y)
  list(data_reefs_sf = data_reefs_sf, data_reefs_df = data_reefs_df)
}

calculate_reef_hcc <- function(spatial_grid, spde, all_effects_hcc, data_reefs_df, data_reefs_sf, reefs_poly_sf) {
  
  spatial_grid_pts_df <- spatial_grid |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    arrange(Longitude, Latitude)

  data_reefs_sample_hcc <- inla.mesh.project(spde$mesh,
    loc = as.matrix(data_reefs_df[, 1:2]),
    all_effects_hcc
  )
  data_reefs_hcc <- data_reefs_sample_hcc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(data_reefs_df) |>
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

  data_reefs_pts_hcc_sf <- data_reefs_hcc |>
    st_as_sf(coords = c("Longitude", "Latitude")) |>
    st_set_crs(st_crs(data_reefs_sf))
  sf_use_s2(FALSE)
  data_reefs_pts_hcc_sf <- data_reefs_pts_hcc_sf |>
    st_intersection(reefs_poly_sf)
  sf_use_s2(TRUE)
  list(data_reefs_sample_hcc = data_reefs_sample_hcc,
    data_reefs_hcc =  data_reefs_hcc,
    data_reefs_pts_hcc_sf = data_reefs_pts_hcc_sf)
}

calculate_reef_sc <- function(spatial_grid, spde, all_effects_sc, data_reefs_df, data_reefs_sf, reefs_poly_sf) {
  
  spatial_grid_pts_df <- spatial_grid |>
    st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y) |>
    arrange(Longitude, Latitude)

  data_reefs_sample_sc <- inla.mesh.project(spde$mesh,
    loc = as.matrix(data_reefs_df[, 1:2]),
    all_effects_sc
  )
  data_reefs_sc <- data_reefs_sample_sc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(data_reefs_df) |>
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

  data_reefs_pts_sc_sf <- data_reefs_sc |>
    st_as_sf(coords = c("Longitude", "Latitude")) |>
    st_set_crs(st_crs(data_reefs_sf))
  sf_use_s2(FALSE)
  data_reefs_pts_sc_sf <- data_reefs_pts_sc_sf |>
    st_intersection(reefs_poly_sf)
  sf_use_s2(TRUE)
  list(data_reefs_sample_sc = data_reefs_sample_sc,
    data_reefs_sc =  data_reefs_sc,
    data_reefs_pts_sc_sf = data_reefs_pts_sc_sf)
}

calculate_reef_ma <- function(data_reefs_hcc, data_reefs_sc, data_reefs_sf, reefs_poly_sf) {
  data_reefs_ma <- data_reefs_hcc |>
    rename(HCC = Value) |>
    full_join(data_reefs_sc |> rename(SC = Value)) |>
    mutate(
      Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
      MA = Total_Avail,
      Value = qlogis(MA)
    ) |>
    dplyr::select(-HCC, -SC, -Total_Avail, -MA)

  data_reefs_pts_ma_sf <- data_reefs_ma |>
    st_as_sf(coords = c("Longitude", "Latitude")) |>
    st_set_crs(st_crs(data_reefs_sf))
  sf_use_s2(FALSE)
  data_reefs_pts_ma_sf <- data_reefs_pts_ma_sf |>
    st_intersection(reefs_poly_sf)
  sf_use_s2(TRUE)

  list(
    data_reefs_ma =  data_reefs_ma,
    data_reefs_pts_ma_sf = data_reefs_pts_ma_sf)
}

combine_reef_benthos <- function(data_reefs_pts_hcc_sf, data_reefs_pts_sc_sf, data_reefs_pts_ma_sf) {
  data_reefs_pts_hcc_sf |>
  rename(HCC = Value) |>
  bind_cols(data_reefs_pts_sc_sf |>
              dplyr::select(SC = Value) |> st_drop_geometry()) |>
  bind_cols(data_reefs_pts_ma_sf |>
              dplyr::select(MA = Value) |> st_drop_geometry())
}
