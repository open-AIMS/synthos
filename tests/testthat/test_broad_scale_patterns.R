test_that("Create SPDE mesh works", {
  library(sf)
  library(INLA)

  config <- list(
    seed = 1,
    crs = 4326,
    alpha = 2,
    kappa = 1
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

  mesh <- create_spde_mesh(spatial_grid, config)
  
  testthat::expect(
    inherits(mesh, c("inla.mesh")),
    "mesh must be an inla.mesh object"
  )
})

