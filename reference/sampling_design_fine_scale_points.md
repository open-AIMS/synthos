# Generate Fine-Scale Photo-Point Observations

Fine-Scale Photo-Point Sampling

## Usage

``` r
sampling_design_fine_scale_points(data_fixed_locs_obs, config_pt)
```

## Arguments

- data_fixed_locs_obs:

  An `sf` object representing fine-scale transect observations produced
  by
  [`sampling_design_fine_scale_fixed()`](https://open-aims.github.io/synthos/reference/sampling_design_fine_scale_fixed.md).
  Must include:

  - `Reef`, `Site`, `Transect`, `Year`, `Date`

  - Benthic cover values: `HCC`, `SC`, `MA`

  - Transect-level coordinates

- config_pt:

  A list containing fine-scale point sampling parameters:

  - `Depths` – number of depths per transect

  - `Depth_effect_multiplier` – magnitude of the depth effect

  - `Number_of_transects_per_site` – number of transects (used for
    checks)

  - `Number_of_frames_per_transect` – frames per transect

  - `Points_per_frame` – points per frame

  - `seed` – optional random seed (if provided)

## Value

An `sf` object representing point-level sampling, containing:

- Reef, Site, Transect, Year, Depth, Frame, Point Number

- Group (HCC, SC, MA)

- Value – benthic cover value after depth effects

- geometry – inherited from transect-level coordinates

## Details

Distributes benthic cover (HCC, SC, MA) across photo-frames and
point-level observations. Expands fine-scale transect data to include
depth effects, allocates points in proportion to benthic cover, and
returns a point-level sampling dataset suitable for analysis or
modelling.

This function:

- Expands each transect across multiple depths

- Applies a randomised depth effect to benthic cover values

- Converts cover percentages to point allocations

- Expands into individual point observations and assigns frame/point IDs

## Author

Murray
