# Combine Reef-level Disturbance Data with Benthos

Combine Reef-level Disturbance Data

## Usage

``` r
combine_reef_disturbances(
  benthos_reefs_pts,
  data_reefs_pts_cyc_sf,
  data_reefs_pts_dhw_sf,
  data_reefs_pts_other_sf
)
```

## Arguments

- benthos_reefs_pts:

  A data frame containing reef-level benthos data (HCC, SC, MA).

- data_reefs_pts_cyc_sf:

  An sf object containing reef-level cyclonic disturbance values.

- data_reefs_pts_dhw_sf:

  An sf object containing reef-level degree heating week (DHW) values.

- data_reefs_pts_other_sf:

  An sf object containing reef-level other disturbance values.

## Value

A data.frame containing reef-level benthos and disturbance values with
columns: `HCC`, `SC`, `MA`, `CYC`, `DHW`, `OTHER`, along with reef
coordinates.

## Details

Combines reef-level disturbance data (CYC, DHW, OTHER) with the benthos
data frame, preserving reef coordinates and removing geometry columns
from disturbance data.

## Author

Murray
