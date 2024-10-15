##' Create SPDE mesh
##'
##' This function creates a SPDE mesh for a given spatial grid and configuration.
##' @title Create SPDE mesh
##' @param spatial_grid 
##' @param config 
##' @return mesh 
##' @author Murray
create_spde_mesh <- function(spatial_grid, config) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  mesh_pars <- c(1, 0.5, 0.1, 1, 0.5) *
    sqrt(config$alpha - ncol(spatial_grid_pts_df) / 2) / config$kappa
  s <- INLA::inla.mesh.segment(
    spatial_grid_pts_df[chull(spatial_grid_pts_df), ]
  )
  mesh <- INLA::inla.mesh.2d(
    spatial_grid_pts_df[chull(spatial_grid_pts_df), ],
    max.edge = mesh_pars[1:2],
    cutoff = mesh_pars[3],
    offset = mesh_pars[4:5],
    boundary = s
  )
  return(mesh)
}

create_spde_matern <- function(spatial_grid, mesh, config) {
  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  spde <- INLA::inla.spde2.matern(mesh, alpha = config$alpha)
  ## calculate precision matrix from the parameter values (theta)
  theta <- c(-0.5 * log(4 * pi * config$variance * config$kappa^2), log(config$kappa))
  Q <- INLA::inla.spde2.precision(spde, theta = theta)
  ## calculate a lattic projection to and from mesh
  A <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(spatial_grid_pts_df)
  )
  # OR
  ## A <- inla.mesh.project(mesh = mesh,
  ##                        loc = as.matrix(spatial_grid_pts_df ))$A
  list(spde = spde, Q = Q, A = A)
}

##' Create a mesh and SPDE
##'
##' This function creates a mesh and SPDE for a given spatial grid and
##' configuration.  It is primarily a wrapper for the create_spde_mesh and
##' create_spde_matern functions.
##' @title Create a mesh and SPDE
##' @param spatial_grid
##' A sfc POINT object representing the full spatial grid/
##' @param config
##' A list that should contains the following parameters:
##' - alpha: the smoothness parameter of the Matern covariance function
##' - kappa: the range parameter of the Matern covariance function
##' - variance: the variance parameter of the Matern covariance function
##' @return a list containing the mesh, spde, Q, and A
##' @examples
##' library(sf)
##' library(INLA)
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
##' @author Murray
##' @export
create_spde <- function(spatial_grid, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("alpha", "kappa", "variance")),
    sort(names(config))
  )
  mesh <- create_spde_mesh(spatial_grid, config)
  spde <- create_spde_matern(spatial_grid, mesh, config)
  list(mesh = mesh, spde = spde$spde, Q = spde$Q, A = spde$A)
}
