# Create Synthetic Reef Landscape with Benthos and Disturbances

Create Synthetic Reef Landscape

## Usage

``` r
create_synthetic_reef_landscape(
  spatial_grid,
  config_sp,
  include_disturbances = FALSE,
  verbose = FALSE
)
```

## Arguments

- spatial_grid:

  An sf object representing the spatial grid over which to simulate
  reefs.

- config_sp:

  A list of config_spuration parameters controlling the simulation,
  including: - `seed`: random seed - `model`, `psill`, `range`,
  `nugget`: variogram/SPDE parameters - `patch_threshold`, `reef_width`:
  reef patch parameters - `years`: numeric vector of years to simulate -
  `dhw_weight`, `cyc_weight`, `other_weight`: weights for disturbance
  layers - `hcc_cover_range`, `hcc_growth`: baseline range and growth
  rate of hard coral - `sc_cover_range`, `sc_growth`: baseline range and
  growth rate of soft coral

- include_disturbances:

  Logical; if TRUE, reef-level disturbances (CYC, DHW, OTHER) are
  included in the output.

- verbose:

  Logical; if TRUE, prints progress messages during simulation.

## Value

A data.frame containing reef-level benthos data (`HCC`, `SC`, `MA`) and,
if `include_disturbances = TRUE`, disturbance values (`CYC`, `DHW`,
`OTHER`), along with reef coordinates.

## Details

Generates a synthetic reef landscape from a spatial grid and
config_spuration parameters. This includes generating synthetic fields,
patches, and reefs, SPDE projections, baseline and synthetic benthos
layers (hard coral, soft coral, macroalgae), and optionally reef-level
disturbance layers (CYC, DHW, OTHER). Reef coordinates are preserved in
the output data frame.

## Author

Murray
