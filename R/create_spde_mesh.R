#' Create an SPDE Mesh
#'
#' Constructs a 2D INLA SPDE mesh from the convex hull of an input spatial grid,
#' using mesh parameters derived from the model configuration.
#'
#' @title Create SPDE Mesh
#'
#' @param spatial_grid An sf object representing the spatial grid from which
#'   mesh nodes will be generated.
#' @param config A list containing SPDE parameters (e.g., `alpha`, `kappa`)
#'   used to define mesh spacing and scaling.
#'
#' @return An `inla.mesh` object representing the constructed SPDE mesh.
#'
#' @author Murray
#' @export
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
    cutoff  = mesh_pars[3],
    offset  = mesh_pars[4:5],
    boundary = s
  )

  return(mesh)
}
