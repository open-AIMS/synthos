# Create SPDE Matern Components

Constructs the SPDE Matern model components (spde object, precision
matrix, and projection matrix) for a given spatial grid and INLA mesh.

## Usage

``` r
create_spde_matern(spatial_grid, mesh, config)
```

## Arguments

- spatial_grid:

  An sf object containing the spatial grid.

- mesh:

  An INLA mesh object.

- config:

  A list with SPDE parameters, including `alpha`, `variance`, and
  `kappa`.

## Value

A list containing:

- spde:

  The INLA SPDE Matern model.

- Q:

  The precision matrix derived from SPDE parameters.

- A:

  The projection matrix mapping mesh nodes to grid points.

## Author

Murray
