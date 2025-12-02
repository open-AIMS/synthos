# Fine-Scale Percent Cover

Fine-Scale Percent Cover Estimates

## Usage

``` r
sampling_design_fine_scale_cover(data_fixed_locs_obs, config_pt)
```

## Arguments

- data_fixed_locs_obs:

  An `sf` object representing fine-scale transect observations generated
  by
  [`sampling_design_fine_scale_fixed()`](https://open-aims.github.io/synthos/reference/sampling_design_fine_scale_fixed.md).
  Must contain:

  - `Reef`, `Site`, `Transect`, `Year`, `Date`

  - `HCC`, `SC`, `MA` – benthic cover values

  - geometry associated with transects

- config_pt_lrge:

  A list containing quadrat-level sampling parameters:

  - `Depths` – number of depths per transect

  - `Depth_effect_multiplier` – magnitude of depth effects

  - `Number_of_quadrats_per_transect` – number of quadrats

  - `Quad_sigma` – standard deviation of quadrat-level noise

  - `seed` – optional random seed

## Value

An `sf` object containing quadrat-level percent cover values with:

- Depth, Quadrats (`Quad`)

- Group (HCC, SC, MA)

- Percent cover values (`Value`) after depth + quadrat effects

- all original transect identifiers and geometry

## Details

Generates fine-scale percent cover estimates by expanding benthic cover
values across quadrats. Applies depth effects, adds quadrat-level
variation, and returns an `sf` object with quadrat-based percent cover
suitable for downstream modelling and analysis.

This function:

- expands transects across depth levels

- reshapes benthic groups into long format

- adds a depth effect using a scaled random component

- generates quadrats and applies quadrat-level noise on the logit scale

- converts back to percent cover (0–100%)

## Author

Murray
