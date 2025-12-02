# Create SPDE Mesh

Create an SPDE Mesh

## Usage

``` r
create_spde_mesh(spatial_grid, config)
```

## Arguments

- spatial_grid:

  An sf object representing the spatial grid from which mesh nodes will
  be generated.

- config:

  A list containing SPDE parameters (e.g., `alpha`, `kappa`) used to
  define mesh spacing and scaling.

## Value

An `inla.mesh` object representing the constructed SPDE mesh.

## Details

Constructs a 2D INLA SPDE mesh from the convex hull of an input spatial
grid, using mesh parameters derived from the model configuration.

## Author

Murray
