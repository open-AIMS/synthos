# Combine Reef-level Benthos Data

Combine Reef-level Benthos Data

## Usage

``` r
combine_reef_benthos(
  data_reefs_pts_hcc_sf,
  data_reefs_pts_sc_sf,
  data_reefs_pts_ma_sf
)
```

## Arguments

- data_reefs_pts_hcc_sf:

  An sf object containing reef-level hard coral cover values.

- data_reefs_pts_sc_sf:

  An sf object containing reef-level soft coral cover values.

- data_reefs_pts_ma_sf:

  An sf object containing reef-level macroalgae cover values.

## Value

A data.frame containing reef-level benthos values with columns `HCC`,
`SC`, `MA`, along with reef coordinates (`Longitude`, `Latitude`) and
any other original metadata.

## Details

Combines reef-level benthos data (hard coral, soft coral, and
macroalgae) into a single data frame, preserving reef coordinates and
removing geometry.

## Author

Murray
