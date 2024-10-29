#devtools::test(filter = "sampling_design_fixed")
test_that("Generate large scale fixed design", {
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
  config <- list(n_locs = 25, n_sites = 2, seed = 123)
  benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(benthos_reefs_pts, config)
  testthat::expect(
    inherits(benthos_fixed_locs_sf, "sf"),
    "benthos_fixed_locs_sf should be a sf"
  )
  testthat::expect_contains(
    names(benthos_fixed_locs_sf),
    c("Reef", "geometry", "Site", "Year", "HCC", "SC", "MA")
  )
})

test_that("Generate fine scale fixed design", {
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
  config <- list(n_locs = 25, n_sites = 2, seed = 123)
  benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(benthos_reefs_pts, config)
  config <- list(
    years =  1:12,
    Number_of_transects_per_site = 5,
    Depths = 2,
    Number_of_frames_per_transect = 100,
    Points_per_frame = 5,
    ## Note, the following are on the link scale
    hcc_site_sigma = 0.5, # variability in Sites within Locations
    hcc_transect_sigma = 0.2, # variability in Transects within Sites
    hcc_sigma = 0.1, # random noise
    sc_site_sigma = 0.05, # variability in Sites within Locations
    sc_transect_sigma = 0.02, # variability in Transects within Sites
    sc_sigma = 0.01, # random noise
    ma_site_sigma = 0.5, # variability in Sites within Locations
    ma_transect_sigma = 0.2, # variability in Transects within Sites
    ma_sigma = 0.1 # random noise
  )
  benthos_fixed_locs_obs <- sampling_design_fine_scale_fixed(benthos_fixed_locs_sf, config)
  testthat::expect(
    inherits(benthos_fixed_locs_obs, "sf"),
    "benthos_fixed_locs_obs should be a sf"
  )
  testthat::expect_contains(
    names(benthos_fixed_locs_obs),
    c("Reef", "geometry", "Site", "Transect", "Year", "HCC", "SC", "MA")
  )
})

test_that("Generate fine scale points for fixed design", {
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
  config <- list(n_locs = 25, n_sites = 2, seed = 123)
  benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(benthos_reefs_pts, config)
  config <- list(
    years =  1:12,
    Number_of_transects_per_site = 5,
    Depths = 2,
    Number_of_frames_per_transect = 100,
    Points_per_frame = 5,
    ## Note, the following are on the link scale
    hcc_site_sigma = 0.5, # variability in Sites within Locations
    hcc_transect_sigma = 0.2, # variability in Transects within Sites
    hcc_sigma = 0.1, # random noise
    sc_site_sigma = 0.05, # variability in Sites within Locations
    sc_transect_sigma = 0.02, # variability in Transects within Sites
    sc_sigma = 0.01, # random noise
    ma_site_sigma = 0.5, # variability in Sites within Locations
    ma_transect_sigma = 0.2, # variability in Transects within Sites
    ma_sigma = 0.1 # random noise
  )
  benthos_fixed_locs_obs <- sampling_design_fine_scale_fixed(benthos_fixed_locs_sf, config)
  config <- list(
    Depths = 2,
    Depth_effect_multiplier = 2,
    Number_of_transects_per_site = 5,
    Number_of_frames_per_transect = 100,
    Points_per_frame = 5
  )
  benthos_fixed_locs_points <- sampling_design_fine_scale_points(benthos_fixed_locs_obs, config)
  testthat::expect(
    inherits(benthos_fixed_locs_points, "sf"),
    "benthos_fixed_locs_points should be a sf"
  )
  testthat::expect_contains(
    names(benthos_fixed_locs_points),
    c("Reef", "geometry", "Site", "Transect", "Year", "HCC", "SC", "MA")
  )
})
