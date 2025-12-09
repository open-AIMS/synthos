#' Generate Simulation Settings for Synthetic Reef Landscapes
#'
#' Creates and saves all configuration lists required to simulate synthetic
#' reef landscapes, including spatial, large-scale, fine-scale, and
#' point-based sampling parameters. Each configuration list is also assigned to
#' the global environment for use by downstream modelling functions.
#'
#' @title Generate Settings for Synthetic Reef Simulations
#'
#' @param nreefs Integer. Number of reef locations to simulate.
#' @param nsites Integer. Number of sites per reef to simulate.
#' @param nyears Integer. Number of years to simulate.
#'
#' @details
#' **Spatio-temporal parameters (config_sp)**
#' \itemize{
#'   \item \strong{seed}: Random seed controlling all stochastic processes.
#'   \item \strong{crs}: Coordinate reference system (EPSG code) applied to
#'     the spatial domain.
#'   \item \strong{model}: Variogram model used to generate spatial
#'     covariance. Supported models include "Sph", "Exp", "Gau", "Lin", "Mat",
#'     "Ste", "Pen", "Hug", "Hol", "Cor", and composite forms such as
#'     "Sphlin", "Sphexp", "Sphgaus", "Sphmat", "Sphste", "Sphpen", "Sphhug",
#'     "Sphhol", "Sphcor".
#'   \item \strong{psill}: Partial sill. The variance explained by spatial
#'     structure; the difference between the sill (asymptotic variance of the
#'     variogram) and the nugget. Represents the plateau reached as lag
#'     distance increases.
#'   \item \strong{range}: Distance at which spatial correlation becomes
#'     negligible. Beyond this range, locations are effectively
#'     uncorrelated.
#'   \item \strong{nugget}: Variance at zero distance (the variogram
#'     intercept). Captures measurement error and micro-scale variability
#'     below the sampling resolution.
#'   \item \strong{alpha}: Smoothness parameter of the spatial field. Larger
#'     values produce smoother fields; smaller values yield more rugged
#'     surfaces. In SPDE models, controls the order of the differential
#'     operator.
#'   \item \strong{kappa}: Controls the spatial scale (range) of the SPDE
#'     spatial field. Smaller values correspond to broader spatial
#'     correlation; larger values produce shorter-range correlation.
#'   \item \strong{variance}: Variance parameter used in the construction of
#'     the Matérn precision matrix for the spatial field.
#'   \item \strong{patch\_threshold}: Numeric cutoff below which the spatial
#'     field is masked, producing discrete habitat patches.
#'   \item \strong{reef\_width}: Half-width of the reef "ribbon" used to
#'     represent benthic habitat around the patch outline.
#'   \item \strong{years}: Sequence of years included in the simulation.
#'   \item \strong{dhw\_weight}: Weight determining the relative influence of
#'     degree heating weeks (thermal stress) on benthic cover.
#'   \item \strong{cyc\_weight}: Weight determining the influence of cyclone
#'     disturbances on benthic cover.
#'   \item \strong{other\_weight}: Weight for additional disturbances not
#'     explicitly modelled.
#'   \item \strong{hcc\_cover\_range}: Minimum and maximum expected hard coral
#'     cover over space and time.
#'   \item \strong{hcc\_growth}: Annual growth rate of hard coral cover.
#'   \item \strong{sc\_cover\_range}: Minimum and maximum expected soft coral
#'     cover over space and time.
#'   \item \strong{sc\_growth}: Annual growth rate of soft coral cover.
#' }
#'
#' **Large-scale sampling parameters (config_lrge)**
#' \itemize{
#'   \item \strong{n_locs}: Number of reef locations.
#'   \item \strong{n_sites}: Number of sites per reef.
#'   \item \strong{seed}: Seed for sampling design reproducibility.
#' }
#'
#' **Fine-scale sampling parameters (config_fine)**  
#' (Site, transect, and random-effect structure for each benthic group)
#' \itemize{
#'   \item \strong{years}: Years included in the fine-scale simulation.
#'   \item \strong{Number_of_transects_per_site}: Number of transects simulated at each site.
#'   \item \strong{Depths}: Number of depth strata simulated.
#'   \item \strong{hcc_site_sigma}, \strong{hcc_transect_sigma},
#'     \strong{hcc_sigma}: Variances for site-level, transect-level, and
#'     residual random effects for hard coral cover.
#'   \item \strong{sc_site_sigma}, \strong{sc_transect_sigma},
#'     \strong{sc_sigma}: Equivalent random-effect variances for soft coral.
#'   \item \strong{ma_site_sigma}, \strong{ma_transect_sigma},
#'     \strong{ma_sigma}: Equivalent random-effect variances for macroalgae.
#' }
#'
#' **Point-based sampling parameters (config_pt)**  
#' Defines photo-quadrat and point-count structure.
#' \itemize{
#'   \item \strong{Depths}: Number of depth strata.
#'   \item \strong{Depth_effect_multiplier}: Multiplier controlling depth-related differences
#'     in benthic cover.
#'   \item \strong{Number_of_transects_per_site}: Transects per site.
#'   \item \strong{Number_of_frames_per_transect}: Number of photo frames per transect.
#'   \item \strong{Number_of_quadrats_per_transect}: Quadrat count per transect.
#'   \item \strong{Points_per_frame}: Number of random points per photographic frame.
#'   \item \strong{Quad_sigma}: Random variation among quadrats.
#' }
#'
#' @return Invisibly returns a list containing the four configuration lists:
#'   \code{config_sp}, \code{config_lrge}, \code{config_fine}, and
#'   \code{config_pt}. Each is also assigned to the global environment.
#'
#' @author Murray
#'
#' @export
generateSettings <- function(nreefs, nsites, nyears){

## Config of the spatio-temporal model
config_sp <- list(
  seed = 1,
  crs = 4326,
  model = "Exp",
  psill = 1,
  range = 15,
  nugget = 0,
  alpha = 2,
  kappa = 1,
  variance = 1,
  patch_threshold = 1.75,
  reef_width = 0.01,
  years = 1:nyears,
  dhw_weight = 0.8,
  cyc_weight = 0.19,
  other_weight = 0.01,
  hcc_cover_range = c(0.1, 0.7),
  hcc_growth = 0.3,
  sc_cover_range = c(0.01, 0.1),
  sc_growth =  0.3
)
assign("config_sp", config_sp, envir = .GlobalEnv)

## Config of sampling design for large scale details
config_lrge <- list(n_locs = nreefs, n_sites = nsites, seed = 123)
assign("config_lrge", config_lrge, envir = .GlobalEnv)

## Config for sampling details for fine scale details 
config_fine <- list(
  years =  1:nyears,
  Number_of_transects_per_site = 5,
  Depths = 2,
  ## Note, the following are on the link scale
  hcc_site_sigma = 0.5, # variability in Sites within Locations
  hcc_transect_sigma = 0.2, # variability in Transects within Sites
  hcc_sigma = 0.1, # random noise

  sc_site_sigma = 0.05, # variability in Sites within Locations
  sc_transect_sigma = 0.02, # variability in Transects within Sites
  sc_sigma = 0.01, # random noise

  ma_site_sigma = 0.5, # variability in Sites within Locations
  ma_transect_sigma = 0.2, # variability in Transects within Sites
  ma_sigma = 0.1 # random noise
)
assign("config_fine", config_fine, envir = .GlobalEnv)

## Generate point-based data 
config_pt <- list(
  Depths = 2,
  Depth_effect_multiplier = 2,
  Number_of_transects_per_site = 5,
  Number_of_frames_per_transect = 100,
  Number_of_quadrats_per_transect = 10,
  Points_per_frame = 50,
  Quad_sigma = 0.5
)
assign("config_pt", config_pt, envir = .GlobalEnv)

config_list <- list(config_sp, config_lrge, config_fine, config_pt)
}
