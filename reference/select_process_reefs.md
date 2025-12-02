# Select and Process Synthetic Reefs

Select Reefs and Prepare Synthetic Benthic Data

## Usage

``` r
select_process_reefs(reef_data_synthetic)
```

## Arguments

- reef_data_synthetic:

  A data frame containing synthetic reef survey data. Must include
  columns such as `site_name`, `survey_depth`, `survey_start_date`,
  `survey_transect_number`, and `point_machine_classification`. If
  `data_type = "cover"`, it should also include `cover`.

- data_type:

  A string specifying the type of data to process:

  - `"points"` – process raw point-based benthic classifications

  - `"cover"` – process aggregated percentage cover values

## Value

A processed data frame containing:

- `reef` – reef identifier

- `site` – standardized site name ("Site X")

- `year` – survey year extracted from `survey_start_date`

- `COUNT` – selected point counts (for points data)

- `COVER` – benthic cover (%) for selected reefs

## Details

Processes a synthetic reef dataset to select a subset of reefs for
survey and prepares benthic cover data either at the point level or
aggregated cover level. Reef and site identifiers are extracted, and
site names are standardized.

The function:

- Randomly selects reefs based on `config_lrge$n_locs`

- Aggregates point-based counts or cover values depending on `data_type`

- Extracts year, reef, and site identifiers from `site_name`

- Standardizes site names from "S1" to "Site 1"

- Assigns NA to reefs not selected for survey

## Author

Julie
