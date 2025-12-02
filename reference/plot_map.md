# Benthic COVER_site Heatmap

Plot Benthic COVER_site Heatmap by Reef and Site

## Usage

``` r
plot_map(df)
```

## Arguments

- df:

  A data frame containing:

  - `year` – survey year (numeric)

  - `COVER_site` – benthic COVER_site percentage (0–100)

  - `reef` – reef identifier

  - `site` – site identifier

## Value

A `ggplot` object visualising benthic COVER_site as a heatmap, faceted
by reef.

## Details

Creates a faceted heatmap showing benthic COVER_site (%) through time
for a selected subset of the data. Each tile represents the benthic
COVER_site at a given site and year. The plot is faceted by reef, with
an automatic colour scale and optional custom year breaks for long time
series.

The function:

- uses a viridis colour scale for benthic COVER_site

- applies dynamic year breaks when the time series exceeds 10 years

- facets results by reef for comparison across locations

## Author

Julie
