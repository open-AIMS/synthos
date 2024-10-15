
test_that("Generate DHW disturbance", {
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
  config <- list(years = 1:12, seed = 123)
  dhw <- disturbance_dhw(spatial_grid, matern_projection, config)

  testthat::expect_match(
    names(dhw),
    "dhw_temporal|dhw_effects|dhw_pts_sample|dhw_pts_effects_df",
  )
  testthat::expect(
    inherits(dhw$dhw_temporal, "data.frame"),
    "dhw$dhw_temporal should be a data.frame"
  )
  testthat::expect_contains(names(dhw$dhw_temporal), c("Year", "Y"))

  testthat::expect(
    inherits(dhw$dhw_effects, "matrix"),
    "dhw$dhw_effects should be a matrix"
  )
  testthat::expect(
    inherits(dhw$dhw_pts_sample, "matrix"),
    "dhw$dhw_pts_sample should be a matrix"
  )

  testthat::expect(
    inherits(dhw$dhw_pts_effects_df, "data.frame"),
    "dhw$dhw_pts_effects_df should be a data.frame"
  )
  testthat::expect_contains(names(dhw$dhw_pts_effects_df), c("Longitude", "Latitude", "Year", "Value"))
})


test_that("Generate cyclone disturbance", {
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
  config <- list(years = 1:12, seed = 123)
  cyc <- disturbance_cyc(spatial_grid, matern_projection, config)

  testthat::expect_match(
    names(cyc),
    "cyc_effects|cyc_effects_df|cyc_pts_sample|cyc_pts_effects",
  )

  testthat::expect(
    inherits(cyc$cyc_effects, "matrix"),
    "cyc$cyc_effects should be a matrix"
  )
  testthat::expect(
    inherits(cyc$cyc_effects_df, "data.frame"),
    "cyc$cyc_effects_df should be a data.frame"
  )
  testthat::expect_contains(
    names(cyc$cyc_effects_df),
    c("Latitude", "Longitude", "Y", "Year", "Value")
  )

  testthat::expect(
    inherits(cyc$cyc_pts_sample, "matrix"),
    "cyc$cyc_pts_sample should be a matrix"
  )

  testthat::expect(
    inherits(cyc$cyc_pts_effects, "data.frame"),
    "cyc$cyc_pts_effects should be a data.frame"
  )
  testthat::expect_contains(names(cyc$cyc_pts_effects), c("Longitude", "Latitude", "Year", "Value"))
})

test_that("Generate other disturbance", {
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
  config <- list(years = 1:12, seed = 123)
  other <- disturbance_other(spatial_grid, matern_projection, config)

  testthat::expect_match(
    names(other),
    "other_effects|other_pts_sample|other_pts_effects",
  )

  testthat::expect(
    inherits(other$other_effects, "matrix"),
    "other$other_effects should be a matrix"
  )

  testthat::expect(
    inherits(other$other_pts_sample, "matrix"),
    "other$other_pts_sample should be a matrix"
  )

  testthat::expect(
    inherits(other$other_pts_effects, "data.frame"),
    "other$other_pts_effects should be a data.frame"
  )
  testthat::expect_contains(names(other$other_pts_effects), c("Longitude", "Latitude", "Year", "Value"))
})
