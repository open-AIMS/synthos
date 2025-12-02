# Plot Time Series of Coral Cover by Group and Depth

Plot Coral Cover Time Series by Group

## Usage

``` r
plot_group(df)
```

## Arguments

- df:

  A data frame containing at least the following variables:

  - `year` – survey year (numeric)

  - `COVER` – coral cover proportion (0–1)

  - `reef` – reef identifier

  - `site` – site identifier

  - `survey_depth` – depth category (character or numeric)

  - `point_machine_classification` – benthic classification label

## Value

A `ggplot` object showing coral cover over time, faceted by reef and
coloured by site.

## Details

Generates a faceted time‐series plot of coral cover (%) for a given
subset of the dataset. The plot shows temporal dynamics across sites,
grouped by machine classification (e.g., HCC, SC, MA) and survey depth.
A dynamic colour palette is used so that site colours adjust
automatically to the number of unique sites present in the data slice.

This function:

- creates a dynamic site-level colour palette using
  [`scales::hue_pal()`](https://scales.r-lib.org/reference/pal_hue.html)

- plots coral cover time series at the site level

- facets results across reefs for visual comparison

- automatically labels panels with depth and classification group

## Author

Julie
