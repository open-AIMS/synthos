# Generate Fine-Scale Sampling Design (random)

Fine-Scale Random Sampling Design

## Usage

``` r
sampling_design_fine_scale_random(data_random_locs_sf, config_fine)
```

## Arguments

- data_random_locs_sf:

  An `sf` object representing the large-scale sampling design,
  containing:

  - `Reef` – unique reef identifier

  - `Site` – unique site identifier

  - `geometry` – coordinates for each site

- config_fine:

  A list containing fine-scale sampling and variance parameters,
  including:

  - `Number_of_transects_per_site` – number of transects per site

  - `Depths` – number of depths (not used directly here)

  - `Number_of_frames_per_transect` – frames per transect

  - `Points_per_frame` – points per frame

  - `hcc_site_sigma`, `hcc_transect_sigma`, `hcc_sigma` – random-effect
    SDs for HCC

  - `sc_site_sigma`, `sc_transect_sigma`, `sc_sigma` – random-effect SDs
    for SC

  - `ma_site_sigma`, `ma_transect_sigma`, `ma_sigma` – random-effect SDs
    for MA

## Value

An `sf` data frame representing fine-scale sampling observations,
including reef, site, transect, coordinates, year, and simulated
percentage cover values (`HCC`, `SC`, `MA`).

## Details

Generates a fine-scale sampling hierarchy (transects, frames, and
points) based on a large-scale random sampling design. Returns an `sf`
object with simulated benthic cover values at the transect level,
incorporating site-level and transect-level random effects for HCC, SC,
and MA.

This function:

- Expands large-scale site locations into multiple transects

- Adds site-level and transect-level random effects for HCC, SC, MA

- Converts simulated link-scale values into percentage cover

- Returns a tidy fine-scale dataset ready for analysis or modelling

## Author

Murray
