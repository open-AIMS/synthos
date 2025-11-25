#' Create SPDE Matern Components
#'
#' Constructs the SPDE Matern model components (spde object, precision matrix,
#' and projection matrix) for a given spatial grid and INLA mesh.
#'
#' @param spatial_grid An sf object containing the spatial grid.
#' @param mesh An INLA mesh object.
#' @param config A list with SPDE parameters, including `alpha`, `variance`,
#'   and `kappa`.
#'
#' @return A list containing:
#'   \item{spde}{The INLA SPDE Matern model.}
#'   \item{Q}{The precision matrix derived from SPDE parameters.}
#'   \item{A}{The projection matrix mapping mesh nodes to grid points.}
#'
#' @author Murray
#' @export
create_spde_matern <- function(spatial_grid, mesh, config) {

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  spde <- INLA::inla.spde2.matern(mesh, alpha = config$alpha)

  ## Calculate precision matrix from parameter values (theta)
  theta <- c(
    -0.5 * log(4 * pi * config$variance * config$kappa^2),
    log(config$kappa)
  )
  Q <- INLA::inla.spde2.precision(spde, theta = theta)

  ## Projection matrix to/from mesh
  A <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = as.matrix(spatial_grid_pts_df)
  )

  # Alternative:
  # A <- INLA::inla.mesh.project(mesh, loc = as.matrix(spatial_grid_pts_df))$A

  list(spde = spde, Q = Q, A = A)
}
