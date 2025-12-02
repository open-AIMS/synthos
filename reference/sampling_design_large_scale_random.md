# Large‐scale Randomised Sampling Design

Large‐scale Randomised Sampling Design

## Usage

``` r
sampling_design_large_scale_random(data_reefs_pts_sf, config_lrge)
```

## Arguments

- data_reefs_pts_sf:

  An `sf` object representing the full field. Must contain:

  - `Year` — survey year (numeric)

  - `Reef` — unique reef identifier

  - `HCC`, `SC`, `MA` — benthic cover values (logit scale)

  - `geometry` — spatial geometry

- config_lrge:

  A list with:

  - `n_locs` — number of reef locations to select

  - `n_sites` — number of sites per selected reef

  - `seed` — random seed for reproducibility

## Value

An `sf` object representing the large-scale randomised sampling design.
Includes selected reefs, sampled sites, valid selected years, and all
associated cover values (still on the logit scale).

## Details

Randomly selects reef locations and a fixed number of sites within each
selected location. Valid survey years are then sampled per reef using
temporal constraints (via
[`sample_years_with_condition()`](https://open-aims.github.io/synthos/reference/sample_years_with_condition.md)).

The function:

- Randomly samples reef locations

- Randomly selects a fixed number of sites within each reef

- Uses
  [`sample_years_with_condition()`](https://open-aims.github.io/synthos/reference/sample_years_with_condition.md)
  to select years per reef

- Filters the dataset to retain only valid reef–year combinations

## Author

Julie
