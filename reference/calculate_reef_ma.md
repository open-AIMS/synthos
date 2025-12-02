# Calculate Reef-level Macroalgae Cover

Calculate Reef-level Macroalgae Cover (MA)

## Usage

``` r
calculate_reef_ma(data_reefs_hcc, data_reefs_sc, data_reefs_sf, reefs_poly_sf)
```

## Arguments

- data_reefs_hcc:

  A data.frame containing reef-level hard coral cover (`Value`).

- data_reefs_sc:

  A data.frame containing reef-level soft coral cover (`Value`).

- data_reefs_sf:

  An sf object of the reef sample points.

- reefs_poly_sf:

  An sf polygon object of reef boundaries.

## Value

A list containing:

- `data_reefs_ma` – Reef-level macroalgae cover in long format with
  `Year` and `Value`.

- `data_reefs_pts_ma_sf` – Reef-level macroalgae cover as an sf object.

## Details

Calculates reef-level macroalgae cover (MA) by combining reef-level hard
coral (HCC) and soft coral (SC) cover. Macroalgae is assumed to occupy
the remaining available space, and the function returns reef-level MA
values both as a table and an sf object.

## Author

Murray
