# Generate Settings for Synthetic Reef Simulations

Generate Simulation Settings for Synthetic Reef Landscapes

## Usage

``` r
generateSettings(nreefs, nsites, nyears, dhw_eff, cyc_eff, other_eff)
```

## Arguments

- nreefs:

  Integer. Number of reef locations to simulate.

- nsites:

  Integer. Number of sites to simulate per reef.

- nyears:

  Integer. Number of years to simulate.

- dhw_eff:

  Numeric. Weight controlling the influence of thermal stress (degree
  heating weeks) on benthic cover.

- cyc_eff:

  Numeric. Weight controlling the influence of cyclone disturbances on
  benthic cover.

- other_eff:

  Numeric. Weight for additional disturbances not explicitly modelled.

## Value

Invisibly returns a list containing the four configuration lists:
`config_sp`, `config_lrge`, `config_fine`, and `config_pt`. Each is also
assigned to the global environment.

## Details

Creates and stores all configuration lists required to simulate
synthetic reef landscapes, including spatio-temporal settings,
large-scale sampling structure, fine-scale sampling hierarchy, and
point-based observation parameters. Each configuration list is also
exported to the global environment for downstream modelling functions.

**Spatio-temporal parameters (`config_sp`)**

Controls the spatial domain, covariance structure, temporal sequence,
and disturbance weights.

- **seed**: Random seed controlling all stochastic processes.

- **crs**: Coordinate reference system (EPSG code).

- **model**: Variogram model (e.g., "Exp", "Sph", "Gau", "Mat", etc.).

- **psill**: Partial sill (variance explained by spatial structure).

- **range**: Distance at which spatial correlation becomes negligible.

- **nugget**: Variance at zero distance (micro-scale variability).

- **alpha**: Smoothness of the spatial field.

- **kappa**: Spatial scale parameter in the SPDE representation.

- **variance**: Variance of the Matérn precision matrix.

- **patch_threshold**: Threshold defining habitat patches.

- **reef_width**: Half-width of the simulated reef ribbon.

- **years**: Sequence of simulated years.

- **dhw_weight**: Weight of thermal stress (degree heating weeks).

- **cyc_weight**: Weight of cyclone disturbances.

- **other_weight**: Weight of additional disturbance processes.

- **hcc_cover_range**: Expected range of hard coral cover.

- **hcc_growth**: Annual growth rate of hard coral cover.

- **sc_cover_range**: Expected range of soft coral cover.

- **sc_growth**: Annual growth rate of soft coral cover.

**Large-scale sampling parameters (`config_lrge`)**

- **n_locs**: Number of reef locations.

- **n_sites**: Number of sites per reef.

- **seed**: Seed for sampling reproducibility.

**Fine-scale sampling parameters (`config_fine`)**

Defines site-level, transect-level, and residual variance components for
each benthic group.

- **years**: Years included in the fine-scale simulation.

- **Number_of_transects_per_site**: Transects simulated at each site.

- **Depths**: Number of depth strata.

- **hcc_site_sigma**, **hcc_transect_sigma**, **hcc_sigma**:
  Random-effect variances for hard coral.

- **sc_site_sigma**, **sc_transect_sigma**, **sc_sigma**: Random-effect
  variances for soft coral.

- **ma_site_sigma**, **ma_transect_sigma**, **ma_sigma**: Random-effect
  variances for macroalgae.

**Point-based sampling parameters (`config_pt`)**

Defines the structure of quadrats, frames, and annotation points.

- **Depths**: Number of depth strata.

- **Depth_effect_multiplier**: Strength of depth-related differences.

- **Number_of_transects_per_site**: Transects per site.

- **Number_of_frames_per_transect**: Frames per transect.

- **Number_of_quadrats_per_transect**: Quadrats per transect.

- **Points_per_frame**: Annotation points per frame.

- **Quad_sigma**: Quadrats-level random variation.

## Author

Murray
