# Plot Benthic Cover Dynamics

Generate Benthic Cover Trajectory or Heatmap Plots

## Usage

``` r
plot_synthos(synthos_data, type = "")
```

## Arguments

- synthos_data:

  A data frame containing benthic survey observations. Required
  variables depend on `data_type`:

  - `survey_depth` – depth category

  - `project_name` – project identifier

  - `site_name` – site name (e.g., "Reef10 Site 1")

  - `survey_start_date` – survey timestamp

  - `point_machine_classification` – benthic class (HCC, SC, MA, etc.)
    *used when `data_type = "points"`*

  - `cover` – numeric benthic cover values *used when
    `data_type = "cover"`*

- type:

  A character string specifying the type of plot to generate:

  - `"trajectories"` – time‐series plots (uses
    [`plot_group()`](https://open-aims.github.io/synthos/reference/plot_group.md))

  - `"heatmaps"` – benthic cover heatmaps (uses
    [`plot_map()`](https://open-aims.github.io/synthos/reference/plot_map.md))

- data_type:

  A character string indicating the structure of input data:

  - `"points"` – raw point-based benthic classifications

  - `"cover"` – pre-aggregated benthic cover values

## Value

A named list of `ggplot` objects, where each plot corresponds to a
unique combination of depth and benthic classification group. The
resulting list can be visualised directly or iterated over using
[`purrr::walk()`](https://purrr.tidyverse.org/reference/map.html).

## Details

Produces a set of plots showing temporal patterns in benthic cover
across survey depths and benthic classification groups. The function
supports two data formats: point-level benthic classifications
(`"points"`) and aggregated benthic cover values (`"cover"`). Depending
on the `type` argument, the output is either a series of time‐series
trajectory plots or benthic‐cover heatmaps.

Based on `data_type`, the function:

- aggregates point classifications into benthic cover proportions; or

- averages numeric benthic cover values

Afterwards, it:

- extracts year, reef, and site identifiers

- splits the dataset by survey depth and classification

- passes each subset to either
  [`plot_group()`](https://open-aims.github.io/synthos/reference/plot_group.md)
  or
  [`plot_map()`](https://open-aims.github.io/synthos/reference/plot_map.md)

- returns all plots as a list

## See also

[`plot_group()`](https://open-aims.github.io/synthos/reference/plot_group.md),
[`plot_map()`](https://open-aims.github.io/synthos/reference/plot_map.md)

## Author

Julie
