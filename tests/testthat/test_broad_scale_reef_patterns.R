
test_that("Generate syntheric reef-level hard coral cover", {
  library(sf)
  library(INLA)
  library(dplyr)
  library(ggplot2)
  config <- list(
    seed = 1,
    crs = 4326,
    model = "Exp",
    psill = 1,
    range = 15,
    nugget = 0,
    patch_threshold = 1.75,
    reef_width = 0.01
  )
  spatial_domain <- st_geometry(
    st_multipoint(
      x = rbind(
        c(0, -10),
        c(3, -10),
        c(10, -20),
        c(1, -21),
        c(2, -16),
        c(0, -10)
      )
    )
  ) |>
    st_set_crs(config$crs) |>
    st_cast("POLYGON")
  set.seed(config$seed)
  spatial_grid <- spatial_domain |>
    st_set_crs(NA) |>
    st_sample(size = 10000, type = "regular") |>
    st_set_crs(config$crs)
  simulated_field <- generate_field(spatial_grid, config)
  simulated_patches <- generate_patches(simulated_field, config)
  simulated_reefs <- generate_reefs(simulated_patches, config)
  config <- list(alpha = 2, kappa = 1, variance = 1)
  matern_projection <- create_spde(spatial_grid, config)
  config <- list(years = 1:12, seed = 123)
  dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
  cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
  other <- disturbance_other(spatial_grid, matern_projection, config)
  config <- list(
    years = 1:12, seed = 123,
    dhw_weight = 0.5,
    cyc_weight = 0.4,
    other_weight = 0.1,
    hcc_growth = 0.3,
    sc_growth =  0.3
  )
  all_disturbance_effects <- disturbance_all(
    spatial_grid,
    dhw_effects = dhw$dhw_effects,
    cyc_effects = cyc$cyc_effects,
    other_effects = other$other_effects,
    matern_projection,
    config) 
  baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
  field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
    baseline_hcc$baseline_sample_hcc, matern_projection)
  reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
  reefs_hcc <- calculate_reef_hcc(
    spatial_grid,
    matern_projection,
    field_hcc$all_effects_hcc,
    reefs$data_reefs_df,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf
  )

  testthat::expect_in(
    names(reefs_hcc),
    c("data_reefs_sample_hcc",
    "data_reefs_hcc",
    "data_reefs_pts_hcc_sf")
  )

  testthat::expect(
    inherits(reefs_hcc$data_reefs_sample_hcc, "matrix"),
    "reefs_hcc$data_reefs_sample_hcc should be a matrix"
  )
  testthat::expect(
    inherits(reefs_hcc$data_reefs_hcc, "data.frame"),
    "reefs_hcc$data_reefs_hcc should be a data.frame"
  )
  testthat::expect_contains(
    names(reefs_hcc$data_reefs_hcc),
    c("Longitude", "Latitude", "Year", "Value")
  )
  testthat::expect(
    inherits(reefs_hcc$data_reefs_pts_hcc_sf, "sf"),
    "reefs_hcc$data_reefs_pts_hcc_sf should be a sf"
  )
  testthat::expect_contains(
    names(reefs_hcc$data_reefs_pts_hcc_sf),
    c("Year", "Value", "Reef")
  )
})

test_that("Generate syntheric reef-level soft coral cover", {
  library(sf)
  library(INLA)
  library(dplyr)
  library(ggplot2)
  config <- list(
    seed = 1,
    crs = 4326,
    model = "Exp",
    psill = 1,
    range = 15,
    nugget = 0,
    patch_threshold = 1.75,
    reef_width = 0.01
  )
  spatial_domain <- st_geometry(
    st_multipoint(
      x = rbind(
        c(0, -10),
        c(3, -10),
        c(10, -20),
        c(1, -21),
        c(2, -16),
        c(0, -10)
      )
    )
  ) |>
    st_set_crs(config$crs) |>
    st_cast("POLYGON")
  set.seed(config$seed)
  spatial_grid <- spatial_domain |>
    st_set_crs(NA) |>
    st_sample(size = 10000, type = "regular") |>
    st_set_crs(config$crs)
  simulated_field <- generate_field(spatial_grid, config)
  simulated_patches <- generate_patches(simulated_field, config)
  simulated_reefs <- generate_reefs(simulated_patches, config)
  config <- list(alpha = 2, kappa = 1, variance = 1)
  matern_projection <- create_spde(spatial_grid, config)
  config <- list(years = 1:12, seed = 123)
  dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
  cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
  other <- disturbance_other(spatial_grid, matern_projection, config)
  config <- list(
    years = 1:12, seed = 123,
    dhw_weight = 0.5,
    cyc_weight = 0.4,
    other_weight = 0.1,
    hcc_growth = 0.3,
    sc_growth =  0.3
  )
  all_disturbance_effects <- disturbance_all(
    spatial_grid,
    dhw_effects = dhw$dhw_effects,
    cyc_effects = cyc$cyc_effects,
    other_effects = other$other_effects,
    matern_projection,
    config) 
  baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
  field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
    baseline_sc$baseline_sample_sc, matern_projection)
  reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
  reefs_sc <- calculate_reef_sc(
    spatial_grid,
    matern_projection,
    field_sc$all_effects_sc,
    reefs$data_reefs_df,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf
  )

  testthat::expect_in(
    names(reefs_sc),
    c("data_reefs_sample_sc",
    "data_reefs_sc",
    "data_reefs_pts_sc_sf")
  )

  testthat::expect(
    inherits(reefs_sc$data_reefs_sample_sc, "matrix"),
    "reefs_sc$data_reefs_sample_sc should be a matrix"
  )
  testthat::expect(
    inherits(reefs_sc$data_reefs_sc, "data.frame"),
    "reefs_sc$data_reefs_sc should be a data.frame"
  )
  testthat::expect_contains(
    names(reefs_sc$data_reefs_sc),
    c("Longitude", "Latitude", "Year", "Value")
  )
  testthat::expect(
    inherits(reefs_sc$data_reefs_pts_sc_sf, "sf"),
    "reefs_sc$data_reefs_pts_sc_sf should be a sf"
  )
  testthat::expect_contains(
    names(reefs_sc$data_reefs_pts_sc_sf),
    c("Year", "Value", "Reef")
  )
})

test_that("Generate syntheric reef-level macroalgae cover", {
 library(sf)
 library(INLA)
 library(dplyr)
 library(ggplot2)
 config <- list(
  seed = 1,
  crs = 4326,
  model = "Exp",
  psill = 1,
  range = 15,
  nugget = 0,
  patch_threshold = 1.75,
  reef_width = 0.01
 )
 spatial_domain <- st_geometry(
   st_multipoint(
     x = rbind(
       c(0, -10),
       c(3, -10),
       c(10, -20),
       c(1, -21),
       c(2, -16),
       c(0, -10)
     )
   )
 ) |>
   st_set_crs(config$crs) |>
   st_cast("POLYGON")
 set.seed(config$seed)
 spatial_grid <- spatial_domain |>
   st_set_crs(NA) |>
   st_sample(size = 10000, type = "regular") |>
   st_set_crs(config$crs)
 simulated_field <- generate_field(spatial_grid, config)
 simulated_patches <- generate_patches(simulated_field, config)
 simulated_reefs <- generate_reefs(simulated_patches, config)
 config <- list(alpha = 2, kappa = 1, variance = 1)
 matern_projection <- create_spde(spatial_grid, config)
 config <- list(years = 1:12, seed = 123)
 dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
 cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
 other <- disturbance_other(spatial_grid, matern_projection, config)
 config <- list(
   years = 1:12, seed = 123,
   dhw_weight = 0.5,
   cyc_weight = 0.4,
   other_weight = 0.1,
   hcc_growth = 0.3,
   sc_growth =  0.3
 )
 all_disturbance_effects <- disturbance_all(
   spatial_grid,
   dhw_effects = dhw$dhw_effects,
   cyc_effects = cyc$cyc_effects,
   other_effects = other$other_effects,
   matern_projection,
   config) 
 baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
 field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
  baseline_hcc$baseline_sample_hcc, matern_projection)
 baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
 field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
  baseline_sc$baseline_sample_sc, matern_projection)
 reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
 reefs_hcc <- calculate_reef_hcc(
   spatial_grid,
   matern_projection,
   field_hcc$all_effects_hcc,
   reefs$data_reefs_df,
   reefs$data_reefs_sf,
   simulated_reefs$simulated_reefs_poly_sf
 )
 reefs_sc <- calculate_reef_sc(
   spatial_grid,
   matern_projection,
   field_sc$all_effects_sc,
   reefs$data_reefs_df,
   reefs$data_reefs_sf,
   simulated_reefs$simulated_reefs_poly_sf
 )
 reefs_ma <- calculate_reef_ma(
   reefs_hcc$data_reefs_hcc,
   reefs_sc$data_reefs_sc,
   reefs$data_reefs_sf,
   simulated_reefs$simulated_reefs_poly_sf
 )

  testthat::expect_in(
    names(reefs_ma),
    c("data_reefs_ma",
    "data_reefs_pts_ma_sf")
  )

  testthat::expect(
    inherits(reefs_ma$data_reefs_ma, "data.frame"),
    "reefs_ma$data_reefs_ma should be a data.frame"
  )
  testthat::expect_contains(
    names(reefs_ma$data_reefs_ma),
    c("Longitude", "Latitude", "Year", "Value")
  )
  testthat::expect(
    inherits(reefs_ma$data_reefs_pts_ma_sf, "sf"),
    "reefs_ma$data_reefs_pts_ma_sf should be a sf"
  )
  testthat::expect_contains(
    names(reefs_ma$data_reefs_pts_ma_sf),
    c("Year", "Value", "Reef")
  )
})

test_that("Combine all synthetic reef-level data", {
 library(sf)
 library(INLA)
 library(dplyr)
 library(ggplot2)
 config <- list(
  seed = 1,
  crs = 4326,
  model = "Exp",
  psill = 1,
  range = 15,
  nugget = 0,
  patch_threshold = 1.75,
  reef_width = 0.01
 )
 spatial_domain <- st_geometry(
   st_multipoint(
     x = rbind(
       c(0, -10),
       c(3, -10),
       c(10, -20),
       c(1, -21),
       c(2, -16),
       c(0, -10)
     )
   )
 ) |>
   st_set_crs(config$crs) |>
   st_cast("POLYGON")
 set.seed(config$seed)
 spatial_grid <- spatial_domain |>
   st_set_crs(NA) |>
   st_sample(size = 10000, type = "regular") |>
   st_set_crs(config$crs)
 simulated_field <- generate_field(spatial_grid, config)
 simulated_patches <- generate_patches(simulated_field, config)
 simulated_reefs <- generate_reefs(simulated_patches, config)
 config <- list(alpha = 2, kappa = 1, variance = 1)
 matern_projection <- create_spde(spatial_grid, config)
 config <- list(years = 1:12, seed = 123)
 dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
 cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
 other <- disturbance_other(spatial_grid, matern_projection, config)
 config <- list(
   years = 1:12, seed = 123,
   dhw_weight = 0.5,
   cyc_weight = 0.4,
   other_weight = 0.1,
   hcc_growth = 0.3,
   sc_growth =  0.3
 )
 all_disturbance_effects <- disturbance_all(
   spatial_grid,
   dhw_effects = dhw$dhw_effects,
   cyc_effects = cyc$cyc_effects,
   other_effects = other$other_effects,
   matern_projection,
   config) 
 baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)
 field_hcc <- synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
  baseline_hcc$baseline_sample_hcc, matern_projection)
 baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)
 field_sc <- synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
  baseline_sc$baseline_sample_sc, matern_projection)
 reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
 reefs_hcc <- calculate_reef_hcc(
   spatial_grid,
   matern_projection,
   field_hcc$all_effects_hcc,
   reefs$data_reefs_df,
   reefs$data_reefs_sf,
   simulated_reefs$simulated_reefs_poly_sf
 )
 reefs_sc <- calculate_reef_sc(
   spatial_grid,
   matern_projection,
   field_sc$all_effects_sc,
   reefs$data_reefs_df,
   reefs$data_reefs_sf,
   simulated_reefs$simulated_reefs_poly_sf
 )
 reefs_ma <- calculate_reef_ma(
   reefs_hcc$data_reefs_hcc,
   reefs_sc$data_reefs_sc,
   reefs$data_reefs_sf,
   simulated_reefs$simulated_reefs_poly_sf
 )
 benthos_reefs_pts <- combine_reef_benthos(
   reefs_hcc$data_reefs_pts_hcc_sf,
   reefs_sc$data_reefs_pts_sc_sf,
   reefs_ma$data_reefs_pts_ma_sf
 )

  testthat::expect(
    inherits(benthos_reefs_pts, "data.frame"),
    "reefs_ma$data_reefs_ma should be a data.frame"
  )
  testthat::expect_contains(
    names(reefs_ma$data_reefs_ma),
    c("Longitude", "Latitude", "Year", "Value")
  )
})
