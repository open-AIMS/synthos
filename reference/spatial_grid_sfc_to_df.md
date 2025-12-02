# Convert sf Geometry to Data Frame

Convert sf Geometry to a Coordinate Data Frame

## Usage

``` r
spatial_grid_sfc_to_df(spatial_grid)
```

## Arguments

- spatial_grid:

  An sf object containing point geometries from which coordinates will
  be extracted.

## Value

A data frame with columns `Longitude` and `Latitude`, ordered
lexicographically.

## Details

Extracts point coordinates from an sf geometry column (`sfc`) and
returns them as a tidy data frame with renamed longitude and latitude
columns.

## Author

Murray
