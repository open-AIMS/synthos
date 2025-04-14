##' Calculate baseline hard coral cover
##'
##' Calculate baseline hard coral cover. The baseline represents the
##' spatial pattern of hard coral cover the year prior to sampling.
##' This spatial pattern is defined as a simple sine wave (applied to
##' the centered latitudes) and rotated slightly and projected onto the
##' SPDE grid.
##' \deqn{Y = clong + sin(clat) + 1.5*clong + clat}
##' Note, these values are on the expected link scale (logit).
##' These values are then rescaled to the nominated cover range (on the
##' logit scale).
##' @title Calculate baseline hard coral cover 
##' @param spatial_grid 
##' A sfc POINT object representing the full spatial grid/
##' @param spde 
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @param cover_range 
##' A numeric vector of length 2 that represents the broad scale range
##' of the coral cover.  The range needs to be between 0 and 1 but cannot
##' include 0 (0%) or 1 (%).  The default values are c(0.1, 0.7) equating
##' to a range of 10% to 70% cover across the spatial domain.
##' @return A list containing the baseline sample, baseline effects,
##' baseline points sample and baseline points effects
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' ggplot(baseline_hcc$baseline_pts_effects_hcc, aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = Value)) +
##'   scale_fill_gradientn(colors = terrain.colors(10)) +
##'   coord_sf(crs = 4236) +
##'   theme_bw(base_size = 7) +
##'   theme(
##'     axis.title = element_blank(),
##'     legend.position = c(0.95, 0.95),
##'     legend.justification = c(1, 1)
##'   )
##' ## Now on a response (%) scale
##' ggplot(baseline_hcc$baseline_pts_effects_hcc, aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = 100 * plogis(Value))) +
##'   scale_fill_gradientn("Cover (%)", colors = terrain.colors(10)) +
##'   coord_sf(crs = 4236) +
##'   theme_bw(base_size = 7) +
##'   theme(
##'     axis.title = element_blank(),
##'     legend.position = c(0.95, 0.95),
##'     legend.justification = c(1, 1)
##'   )
##' @export
baseline_hard_coral_cover <- function(spatial_grid, spde, cover_range = c(0.1, 0.7)) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  testthat::expect(
    min(cover_range) > 0 & max(cover_range) < 1,
    "cover range must be between 0 and 1"
  )
  cover_range <- qlogis(cover_range)

  baseline_sample_hcc <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    dplyr::mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = clong + sin(clat) + # rnorm(1,0,1) +
        1.5 * clong + clat
    ) |>
    dplyr::mutate(Y = scales::rescale(Y, to = cover_range))

  baseline_effects_hcc <- baseline_sample_hcc |>
    dplyr::select(Y) |>
    as.matrix()
  baseline_pts_sample_hcc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    baseline_effects_hcc
  )
  baseline_pts_effects_hcc <- baseline_pts_sample_hcc |>
    cbind() |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = as.numeric(Year))

  list(baseline_sample_hcc = baseline_sample_hcc,
    baseline_effects_hcc = baseline_effects_hcc,
    baseline_pts_sample_hcc = baseline_pts_sample_hcc,
    baseline_pts_effects_hcc = baseline_pts_effects_hcc)
}

##' Create a synthetic field of hard coral cover
##'
##' Create a broad scale synthetic field of hard coral cover. This
##' field is created by adding the baseline field to the effects of
##' the year. The effects of the year are assumed to be on the link
##' scale (logit) and are projected onto the SPDE grid.
##' @title Create a synthetic field of hard coral cover
##' @param spatial_grid
##' A sfc POINT object representing the full spatial grid
##' @param all_effects_df 
##' A data.frame containing the effects of the year
##' @param baseline_sample_hcc 
##' A matrix containing the baseline sample
##' @param spde 
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @return 
##' A list containing the synthetic field, the projected synthetic field
##' and the projected synthetic field as a data.frame
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection, cover_range = c(0.1, 0.7))
##' field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' 
##' field_hcc$all_pts_effects_hcc |>
##'   mutate(Value = plogis(Value)) |> 
##'   group_by(Year) |>
##'   summarise(Mean = mean(Value),
##'     Median = median(Value)) |>
##'   ggplot(aes(x = Year)) +
##'   geom_line(aes(y = Mean), color = "blue") +
##'   geom_line(aes(y = Median), color = "red")   
##' ## Now on a percentage scale
##' field_hcc$all_pts_effects_hcc |>
##'   mutate(Value = plogis(Value)) |> 
##'   ggplot(aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = Value)) +
##'   scale_fill_gradientn(colors = terrain.colors(10)) + 
##'   coord_sf(crs = 4236) +
##'   facet_wrap(~Year, nrow = 2) + 
##'   theme_bw(base_size = 7) +
##'   theme(axis.title = element_blank(),
##'     legend.position = c(0.95,0.95),
##'     legend.justification = c(1,1))
##' @export
synthetic_field_hcc <- function(spatial_grid, all_effects_df, baseline_sample_hcc, spde) {
  testthat::expect(
    inherits(all_effects_df, c("data.frame")),
    "all_effects_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(all_effects_df),
    c("Longitude", "Latitude", "Year", "Y", "Growth_HCC", "Y_HCC")
  )
  testthat::expect(
    inherits(baseline_sample_hcc, c("data.frame")),
    "baseline_sample_hcc must be a data.frame object"
  )
  testthat::expect_contains(
    names(baseline_sample_hcc),
    c("Longitude", "Latitude", "Y")
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  ## Do all this on the link scale so that can use cumsum
  all_effects_hcc <- all_effects_df |>
    dplyr::full_join(baseline_sample_hcc |>
                dplyr::select(Longitude, Latitude, BASE_HCC = Y)) |>
    dplyr::group_by(Longitude, Latitude) |>
    dplyr::mutate(HCC = BASE_HCC + Y_HCC) |>
    dplyr::ungroup() |>
    dplyr::select(-BASE_HCC, -Y_HCC) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = HCC
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  all_pts_sample_hcc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects_hcc
  )
  all_pts_effects_hcc <- all_pts_sample_hcc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Year = config$years[as.numeric(Year)],
      Value = Value
    )
  list(all_effects_hcc = all_effects_hcc,
    all_pts_sample_hcc = all_pts_sample_hcc,
    all_pts_effects_hcc = all_pts_effects_hcc)
}


##' Calculate baseline soft cover
##'
##' Calculate baseline soft cover. The baseline represents the
##' spatial pattern of soft cover the year prior to sampling.
##' This spatial pattern is defined as a simple sine wave (applied to
##' the centered latitudes) and rotated slightly and projected onto the
##' SPDE grid.
##' \deqn{Y = clong + sin(clat) + 1.5*clong - 1.5*clat}
##' Note, these values are on the expected link scale (logit).
##' @title Calculate baseline soft cover 
##' @param spatial_grid 
##' A sfc POINT object representing the full spatial grid/
##' @param spde 
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @param cover_range 
##' A numeric vector of length 2 that represents the broad scale range
##' of the coral cover.  The range needs to be between 0 and 1 but cannot
##' include 0 (0%) or 1 (%).  The default values are c(0.01, 0.1) equating
##' to a range of 0.1% to 10% cover across the spatial domain.
##' @return A list containing the baseline sample, baseline effects,
##' baseline points sample and baseline points effects
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' ggplot(baseline_sc$baseline_pts_effects_sc, aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = Value)) +
##'   scale_fill_gradientn(colors = terrain.colors(10)) +
##'   coord_sf(crs = 4236) +
##'   theme_bw(base_size = 7) +
##'   theme(
##'     axis.title = element_blank(),
##'     legend.position = c(0.95, 0.95),
##'     legend.justification = c(1, 1)
##'   )
##' ## Now on a response (%) scale
##' ggplot(baseline_sc$baseline_pts_effects_sc, aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = 100 * plogis(Value))) +
##'   scale_fill_gradientn("Cover (%)", colors = terrain.colors(10)) +
##'   coord_sf(crs = 4236) +
##'   theme_bw(base_size = 7) +
##'   theme(
##'     axis.title = element_blank(),
##'     legend.position = c(0.95, 0.95),
##'     legend.justification = c(1, 1)
##'   )
##' @export
baseline_soft_coral_cover <- function(spatial_grid, spde, cover_range = c(0.01, 0.1)) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  testthat::expect(
    min(cover_range) > 0 & max(cover_range) < 1,
    "cover range must be between 0 and 1"
  )
  cover_range <- qlogis(cover_range)

  baseline_sample_sc <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::select(Longitude = V1, Latitude = V2) |>
    dplyr::mutate(
      clong = as.vector(scale(Longitude, scale = FALSE)),
      clat = as.vector(scale(Latitude, scale = FALSE)),
      Y = clong + sin(clat) + # rnorm(1,0,1) +
        1.5 * clong + -1.5 * clat
    ) |>
    dplyr::mutate(Y = scales::rescale(Y, to = cover_range))

  baseline_effects_sc <- baseline_sample_sc |>
    dplyr::select(Y) |>
    as.matrix()
  baseline_pts_sample_sc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    baseline_effects_sc
  )
  baseline_pts_effects_sc <- baseline_pts_sample_sc |>
    cbind() |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = config$years[as.numeric(Year)])

  list(baseline_sample_sc = baseline_sample_sc,
    baseline_effects_sc = baseline_effects_sc,
    baseline_pts_sample_sc = baseline_pts_sample_sc,
    baseline_pts_effects_sc = baseline_pts_effects_sc)
}


##' Create a synthetic field of soft coral cover
##'
##' Create a broad scale synthetic field of soft coral cover. This
##' field is created by adding the baseline field to the effects of
##' the year. The effects of the year are assumed to be on the link
##' scale (logit) and are projected onto the SPDE grid.
##' @title Create a synthetic field of soft coral cover
##' @param spatial_grid 
##' A sfc POINT object representing the full spatial grid/
##' @param all_effects_df 
##' A data.frame containing the effects of the year
##' @param baseline_sample_sc 
##' A matrix containing the baseline sample
##' @param spde 
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @return 
##' A list containing the synthetic field, the projected synthetic field
##' and the projected synthetic field as a data.frame
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
##'  baseline_sc$baseline_sample_sc, matern_projection)
##' 
##' field_sc$all_pts_effects_sc |>
##'   mutate(Value = plogis(Value)) |> 
##'   group_by(Year) |>
##'   summarise(Mean = mean(Value),
##'     Median = median(Value)) |>
##'   ggplot(aes(x = Year)) +
##'   geom_line(aes(y = Mean), color = "blue") +
##'   geom_line(aes(y = Median), color = "red")   
##' ## Now on a percentage scale
##' field_sc$all_pts_effects_sc |>
##'   mutate(Value = plogis(Value)) |> 
##'   ggplot(aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = Value)) +
##'   scale_fill_gradientn(colors = terrain.colors(10)) + 
##'   coord_sf(crs = 4236) +
##'   facet_wrap(~Year, nrow = 2) + 
##'   theme_bw(base_size = 7) +
##'   theme(axis.title = element_blank(),
##'     legend.position = c(0.95,0.95),
##'     legend.justification = c(1,1))
##' @export
synthetic_field_sc <- function(spatial_grid, all_effects_df, baseline_sample_sc, spde) {
  testthat::expect(
    inherits(all_effects_df, c("data.frame")),
    "all_effects_df must be a data.frame object"
  )
  testthat::expect_contains(
    names(all_effects_df),
    c("Longitude", "Latitude", "Year", "Y", "Growth_SC", "Y_SC")
  )
  testthat::expect(
    inherits(baseline_sample_sc, c("data.frame")),
    "baseline_sample_sc must be a data.frame object"
  )
  testthat::expect_contains(
    names(baseline_sample_sc),
    c("Longitude", "Latitude", "Y")
  )
  testthat::expect(
    inherits(spde$mesh, c("inla.mesh")),
    "spde$mesh must be a inla.mesh object"
  )
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  ## Do all this on the link scale so that can use cumsum
  all_effects_sc <- all_effects_df |>
    dplyr::full_join(baseline_sample_sc |>
                dplyr::select(Longitude, Latitude, BASE_SC = Y)) |>
    dplyr::group_by(Longitude, Latitude) |>
    dplyr::mutate(SC = BASE_SC + Y_SC) |>
    dplyr::ungroup() |>
    dplyr::select(-BASE_SC, -Y_SC) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = SC
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  all_pts_sample_sc <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects_sc
  )
  all_pts_effects_sc <- all_pts_sample_sc |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(
      Year = config$years[as.numeric(Year)],
      Value = Value
    )
  list(all_effects_sc = all_effects_sc,
    all_pts_sample_sc = all_pts_sample_sc,
    all_pts_effects_sc = all_pts_effects_sc)
}


##' Create a synthetic field of macroalgae cover
##'
##' Create a broad scale synthetic field of macroalgae cover. This
##' field is created by adding the baseline field to the effects of
##' the year.  This is different to hard and soft coral cover in that
##' it is essentially the difference between the total available area
##' and the hard and soft coral - since macroalgae often just expands to
##' fill the space.
##' @param all_effects_df 
##' A data.frame containing the effects of the year
##' @param baseline_sample_ma 
##' A matrix containing the baseline sample
##' @param spde 
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @return 
##' A list containing the synthetic field, the projected synthetic field
##' and the projected synthetic field as a data.frame
##' @author Murray
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
##' field_hcc <- synthetic_field_hcc(all_disturbance_effects$all_effects_df,
##'  baseline_hcc$baseline_sample_hcc, matern_projection)
##' baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
##' field_sc <- synthetic_field_sc(all_disturbance_effects$all_effects_df,
##'  baseline_sc$baseline_sample_sc, matern_projection)
##' field_ma <- synthetic_field_ma(
##'  field_hcc$all_pts_effects_hcc,
##'  field_sc$all_pts_effects_sc
##' )
##' field_ma$all_pts_effects_ma |>
##'   mutate(Value = plogis(Value)) |> 
##'   group_by(Year) |>
##'   summarise(Mean = mean(Value),
##'     Median = median(Value)) |>
##'   ggplot(aes(x = Year)) +
##'   geom_line(aes(y = Mean), color = "blue") +
##'   geom_line(aes(y = Median), color = "red")
##' field_ma$all_pts_effects_ma |>
##'   mutate(Value = plogis(Value)) |> 
##'   ggplot(aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = Value)) +
##'   scale_fill_gradientn(colors = terrain.colors(10)) + 
##'   coord_sf(crs = 4236) +
##'   facet_wrap(~Year, nrow = 2) + 
##'   theme_bw(base_size = 7) +
##'   theme(axis.title = element_blank(),
##'     legend.position = c(0.95,0.95),
##'     legend.justification = c(1,1))
##' @export
synthetic_field_ma <- function(all_pts_effects_hcc, all_pts_effects_sc) {
  ## Do all this on the link scale so that can use cumsum
  all_pts_effects_ma <- all_pts_effects_hcc |>
    dplyr::rename(HCC=Value) |> 
    dplyr::bind_cols(all_pts_effects_sc |>
                dplyr::select(SC=Value)) |>
    dplyr::mutate(Total_Avail = 0.8 - plogis(HCC) + plogis(SC),
      ## MA = Total_Avail*rbeta(n(), 2, 1),
      MA = Total_Avail,
      Value = qlogis(MA)) |>
    dplyr::select(-HCC, -SC, -Total_Avail, -MA)
  list(all_pts_effects_ma = all_pts_effects_ma)
}
