
test_that("Generate baseline hard coral cover", {
  library(sf)
  library(INLA)
  library(ggplot2)
  config <- list(crs = 4326, seed = 123)
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
  config <- list(alpha = 2, kappa = 1, variance = 1)
  matern_projection <- create_spde(spatial_grid, config)

  baseline_hcc <- baseline_hard_coral_cover(spatial_grid, matern_projection)

  testthat::expect_in(
    names(baseline_hcc),
    c("baseline_sample_hcc", "baseline_effects_hcc",
      "baseline_pts_sample_hcc", "baseline_pts_effects_hcc")
  )
  testthat::expect(
    inherits(baseline_hcc$baseline_sample_hcc, "data.frame"),
    "baseline_hcc$baseline_sample_hcc should be a data.frame"
  )
  testthat::expect_contains(
    names(baseline_hcc$baseline_sample_hcc),
    c("Longitude", "Latitude", "Y")
  )
  testthat::expect(
    inherits(baseline_hcc$baseline_effects_hcc, "matrix"),
    "baseline_hcc$baseline_effects_hcc should be a matrix"
  )
  testthat::expect(
    inherits(baseline_hcc$baseline_pts_sample_hcc, "matrix"),
    "baseline_hcc$baseline_pts_sample_hcc should be a matrix"
  )
  testthat::expect(
    inherits(baseline_hcc$baseline_pts_effects_hcc, "data.frame"),
    "baseline_hcc$baseline_pts_effects_hcc should be a data.frame"
  )
  testthat::expect_contains(
    names(baseline_hcc$baseline_pts_effects_hcc),
    c("Longitude", "Latitude", "Year", "Value")
  )
})

test_that("Generate syntheric field hard coral cover", {
  library(sf)
  library(INLA)
  library(dplyr)
  library(ggplot2)
  config <- list(crs=4326, seed = 123)
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

  testthat::expect_in(
    names(field_hcc),
    c("all_effects_hcc", "all_pts_sample_hcc", "all_pts_effects_hcc")
  )
  testthat::expect(
    inherits(field_hcc$all_effects_hcc, "matrix"),
    "field_hcc$all_effects_hcc should be a matrix"
  )
  testthat::expect(
    inherits(field_hcc$all_pts_sample_hcc, "matrix"),
    "field_hcc$all_pts_sample_hcc should be a matrix"
  )
  testthat::expect(
    inherits(field_hcc$all_pts_effects_hcc, "data.frame"),
    "field_hcc$all_pts_effects_hcc should be a data.frame"
  )
  testthat::expect_contains(
    names(field_hcc$all_pts_effects_hcc),
    c("Longitude", "Latitude", "Year", "Value")
  )
})

test_that("Generate baseline soft coral cover", {
  library(sf)
  library(INLA)
  library(ggplot2)
  config <- list(crs = 4326, seed = 123)
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
  config <- list(alpha = 2, kappa = 1, variance = 1)
  matern_projection <- create_spde(spatial_grid, config)

  baseline_sc <- baseline_soft_coral_cover(spatial_grid, matern_projection)

  testthat::expect_in(
    names(baseline_sc),
    c("baseline_sample_sc", "baseline_effects_sc",
      "baseline_pts_sample_sc", "baseline_pts_effects_sc")
  )
  testthat::expect(
    inherits(baseline_sc$baseline_sample_sc, "data.frame"),
    "baseline_sc$baseline_sample_sc should be a data.frame"
  )
  testthat::expect_contains(
    names(baseline_sc$baseline_sample_sc),
    c("Longitude", "Latitude", "Y")
  )
  testthat::expect(
    inherits(baseline_sc$baseline_effects_sc, "matrix"),
    "baseline_sc$baseline_effects_sc should be a matrix"
  )
  testthat::expect(
    inherits(baseline_sc$baseline_pts_sample_sc, "matrix"),
    "baseline_sc$baseline_pts_sample_sc should be a matrix"
  )
  testthat::expect(
    inherits(baseline_sc$baseline_pts_effects_sc, "data.frame"),
    "baseline_sc$baseline_pts_effects_sc should be a data.frame"
  )
  testthat::expect_contains(
    names(baseline_sc$baseline_pts_effects_sc),
    c("Longitude", "Latitude", "Year", "Value")
  )
})

test_that("Generate syntheric field soft coral cover", {
  library(sf)
  library(INLA)
  library(dplyr)
  library(ggplot2)
  config <- list(crs=4326, seed = 123)
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
    sc_growth = 0.3,
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

  testthat::expect_in(
    names(field_sc),
    c("all_effects_sc", "all_pts_sample_sc", "all_pts_effects_sc")
  )
  testthat::expect(
    inherits(field_sc$all_effects_sc, "matrix"),
    "field_sc$all_effects_sc should be a matrix"
  )
  testthat::expect(
    inherits(field_sc$all_pts_sample_sc, "matrix"),
    "field_sc$all_pts_sample_sc should be a matrix"
  )
  testthat::expect(
    inherits(field_sc$all_pts_effects_sc, "data.frame"),
    "field_sc$all_pts_effects_sc should be a data.frame"
  )
  testthat::expect_contains(
    names(field_sc$all_pts_effects_sc),
    c("Longitude", "Latitude", "Year", "Value")
  )
})

test_that("Generate syntheric field macroalgae cover", {
  library(sf)
  library(INLA)
  library(dplyr)
  library(ggplot2)
  config <- list(crs=4326, seed = 123)
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
    sc_growth = 0.3,
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

  field_ma <- synthetic_field_ma(
    field_hcc$all_pts_effects_hcc,
    field_sc$all_pts_effects_sc
  )
  testthat::expect_in(
    names(field_ma),
    c("all_pts_effects_ma")
  )
  testthat::expect(
    inherits(field_ma$all_pts_effects_ma, "data.frame"),
    "field_ma should be a data.frame"
  )
  testthat::expect_contains(
    names(field_ma$all_pts_effects_ma),
    c("Longitude", "Latitude", "Year", "Value")
  )
})

test_that("Create synthetic reef landscale wrapper", {
  library(sf)
  library(stars)
  library(gstat)
  library(INLA)
  config <- list(
    seed = 1,
    crs = 4326,
    model = "Exp",
    psill = 1,
    range = 15,
    nugget = 0,
    alpha = 2,
    kappa = 1,
    variance = 1,
    patch_threshold = 1.75,
    reef_width = 0.01,
    years = 1:12,
    dhw_weight = 0.5,
    cyc_weight = 0.4,
    other_weight = 0.1,
    hcc_growth = 0.3,
    sc_growth =  0.3
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
  benthos_reefs_pts <- create_synthetic_reef_landscape(spatial_grid, config)

  testthat::expect(
    inherits(benthos_reefs_pts, "sf"),
    "benthos_reefs_pts should be a sf"
  )
  testthat::expect_contains(
    names(benthos_reefs_pts),
    c("Year", "HCC", "SC", "MA", "Reef", "geometry")
  )
})
