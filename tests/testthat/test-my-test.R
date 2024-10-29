test_that("Generate simulated field works", {
  library(sf)
  library(gstat)
  config <- list(
    seed = 1,
    crs = 4326,
    model = "Exp",
    psill = 1,
    range = 15,
    nugget = 0
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
  simulated_field <- generate_field(spatial_grid, config) |> suppressMessages() |> suppressWarnings()

  testthat::expect(
    inherits(simulated_field, c("sf", "sfc")),
    "spatial_grid must be an sf object"
  )
  testthat::expect_match(
    names(simulated_field),
    "sim[0-9]*|geometry"
  )
})

test_that("Generate simulated patch works", {
  library(sf)
  library(gstat)
  config <- list(
    seed = 1,
    crs = 4326,
    model = "Exp",
    psill = 1,
    range = 15,
    nugget = 0,
    patch_threshold = 1.75
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
  simulated_field <- generate_field(spatial_grid, config) |> suppressMessages() |> suppressWarnings()
  simulated_patch <- generate_patches(simulated_field, config) |> suppressMessages() |> suppressWarnings()
  testthat::expect(
    inherits(simulated_patch, c("sfc")),
    "simulated_field must be an sfc object"
  )
})

test_that("Generate simulated reefs works", {
  library(sf)
  library(gstat)
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
  simulated_field <- generate_field(spatial_grid, config) |> suppressMessages() |> suppressWarnings()
  simulated_patch <- generate_patches(simulated_field, config) |> suppressMessages() |> suppressWarnings()
  simulated_reefs <- generate_reefs(simulated_patch, config) |> suppressMessages() |> suppressWarnings()
  testthat::expect(
    inherits(simulated_reefs$simulated_reefs_sf, c("sfc")),
    "simulated_field must be an sf object"
  )
  testthat::expect(
    inherits(simulated_reefs$simulated_reefs_poly_sf, c("sf")),
    "simulated_field must be an sf object"
  )
  testthat::expect_in(
    sort(c("geometry", "Reef")),
    sort(colnames(simulated_reefs$simulated_reefs_poly_sf))
  )
})

test_that("Pointify polygons", {
 library(sf)
 library(gstat)
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
 reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)

  testthat::expect_in(
    names(reefs),
    c("data_reefs_sf", "data_reefs_df")
  )
  testthat::expect(
    inherits(reefs$data_reefs_sf, "sf"),
    "reefs$data_reefs_sf should be an sf object"
  )
  testthat::expect(
    inherits(reefs$data_reefs_df, "data.frame"),
    "reefs$data_reefs_df should be an data.frame object"
  )
  testthat::expect_contains(
    names(reefs$data_reefs_df),
    c("Longitude", "Latitude")
  )
})
