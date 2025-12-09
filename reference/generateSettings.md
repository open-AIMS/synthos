# Generate Settings for Synthetic Reef Simulations

Generate Simulation Settings for Synthetic Reef Landscapes

## Usage

``` r
generateSettings(nreefs, nsites, nyears)
```

## Arguments

- nreefs:

  Integer. Number of reef locations to simulate.

- nsites:

  Integer. Number of sites per reef to simulate.

- nyears:

  Integer. Number of years to simulate.

## Value

Invisibly returns a list containing the four configuration lists:
`config_sp`, `config_lrge`, `config_fine`, and `config_pt`. Each is also
assigned to the global environment.

## Details

Creates and saves all configuration lists required to simulate synthetic
reef landscapes, including spatial, large-scale, fine-scale, and
point-based sampling parameters. Each configuration list is also
assigned to the global environment for use by downstream modelling
functions.

**Spatio-temporal parameters (config_sp)**

- **seed**: Random seed controlling all stochastic processes.

- **crs**: Coordinate reference system (EPSG code) applied to the
  spatial domain.

- **model**: Variogram model used to generate spatial covariance.
  Supported models include "Sph", "Exp", "Gau", "Lin", "Mat", "Ste",
  "Pen", "Hug", "Hol", "Cor", and composite forms such as "Sphlin",
  "Sphexp", "Sphgaus", "Sphmat", "Sphste", "Sphpen", "Sphhug", "Sphhol",
  "Sphcor".

- **psill**: Partial sill. The variance explained by spatial structure;
  the difference between the sill (asymptotic variance of the variogram)
  and the nugget. Represents the plateau reached as lag distance
  increases.

- **range**: Distance at which spatial correlation becomes negligible.
  Beyond this range, locations are effectively uncorrelated.

- **nugget**: Variance at zero distance (the variogram intercept).
  Captures measurement error and micro-scale variability below the
  sampling resolution.

- **alpha**: Smoothness parameter of the spatial field. Larger values
  produce smoother fields; smaller values yield more rugged surfaces. In
  SPDE models, controls the order of the differential operator.

- **kappa**: Controls the spatial scale (range) of the SPDE spatial
  field. Smaller values correspond to broader spatial correlation;
  larger values produce shorter-range correlation.

- **variance**: Variance parameter used in the construction of the
  Matérn precision matrix for the spatial field.

- **patch\\threshold**: Numeric cutoff below which the spatial field is
  masked, producing discrete habitat patches.

- **reef\\width**: Half-width of the reef "ribbon" used to represent
  benthic habitat around the patch outline.

- **years**: Sequence of years included in the simulation.

- **dhw\\weight**: Weight determining the relative influence of degree
  heating weeks (thermal stress) on benthic cover.

- **cyc\\weight**: Weight determining the influence of cyclone
  disturbances on benthic cover.

- **other\\weight**: Weight for additional disturbances not explicitly
  modelled.

- **hcc\\cover\\range**: Minimum and maximum expected hard coral cover
  over space and time.

- **hcc\\growth**: Annual growth rate of hard coral cover.

- **sc\\cover\\range**: Minimum and maximum expected soft coral cover
  over space and time.

- **sc\\growth**: Annual growth rate of soft coral cover.

**Large-scale sampling parameters (config_lrge)**

- **n_locs**: Number of reef locations.

- **n_sites**: Number of sites per reef.

- **seed**: Seed for sampling design reproducibility.

**Fine-scale sampling parameters (config_fine)** (Site, transect, and
random-effect structure for each benthic group)

- **years**: Years included in the fine-scale simulation.

- **Number_of_transects_per_site**: Number of transects simulated at
  each site.

- **Depths**: Number of depth strata simulated.

- **hcc_site_sigma**, **hcc_transect_sigma**, **hcc_sigma**: Variances
  for site-level, transect-level, and residual random effects for hard
  coral cover.

- **sc_site_sigma**, **sc_transect_sigma**, **sc_sigma**: Equivalent
  random-effect variances for soft coral.

- **ma_site_sigma**, **ma_transect_sigma**, **ma_sigma**: Equivalent
  random-effect variances for macroalgae.

**Point-based sampling parameters (config_pt)** Defines photo-quadrat
and point-count structure.

- **Depths**: Number of depth strata.

- **Depth_effect_multiplier**: Multiplier controlling depth-related
  differences in benthic cover.

- **Number_of_transects_per_site**: Transects per site.

- **Number_of_frames_per_transect**: Number of photo frames per
  transect.

- **Number_of_quadrats_per_transect**: Quadrat count per transect.

- **Points_per_frame**: Number of random points per photographic frame.

- **Quad_sigma**: Random variation among quadrats.

## Author

Murray
