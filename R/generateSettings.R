#' Generate Simulation Settings for Synthetic Reef Landscapes
#'
#' Creates and stores all configuration lists required to simulate synthetic
#' reef landscapes, including spatio-temporal settings, large-scale sampling
#' structure, fine-scale sampling hierarchy, and point-based observation
#' parameters. Each configuration list is also exported to the global
#' environment for downstream modelling functions.
#'
#' @title Generate Settings for Synthetic Reef Simulations
#'
#' @param nreefs Integer. Number of reef locations to simulate.
#' @param nsites Integer. Number of sites to simulate per reef.
#' @param nyears Integer. Number of years to simulate.
#' @param dhw_eff Numeric. Weight controlling the influence of thermal stress
#'   (degree heating weeks) on benthic cover.
#' @param cyc_eff Numeric. Weight controlling the influence of cyclone 
#'   disturbances on benthic cover.
#' @param other_eff Numeric. Weight for additional disturbances not explicitly
#'   modelled.
#'
#' @details
#' **Spatio-temporal parameters (`config_sp`)**
#'
#' Controls the spatial domain, covariance structure, temporal sequence,
#' and disturbance weights.
#'
#' \itemize{
#'   \item \strong{seed}: Random seed controlling all stochastic processes.
#'   \item \strong{crs}: Coordinate reference system (EPSG code).
#'   \item \strong{model}: Variogram model (e.g., "Exp", "Sph", "Gau", "Mat", etc.).
#'   \item \strong{psill}: Partial sill (variance explained by spatial structure).
#'   \item \strong{range}: Distance at which spatial correlation becomes negligible.
#'   \item \strong{nugget}: Variance at zero distance (micro-scale variability).
#'   \item \strong{alpha}: Smoothness of the spatial field.
#'   \item \strong{kappa}: Spatial scale parameter in the SPDE representation.
#'   \item \strong{variance}: Variance of the Matérn precision matrix.
#'   \item \strong{patch_threshold}: Threshold defining habitat patches.
#'   \item \strong{reef_width}: Half-width of the simulated reef ribbon.
#'   \item \strong{years}: Sequence of simulated years.
#'   \item \strong{dhw_weight}: Weight of thermal stress (degree heating weeks).
#'   \item \strong{cyc_weight}: Weight of cyclone disturbances.
#'   \item \strong{other_weight}: Weight of additional disturbance processes.
#'   \item \strong{hcc_cover_range}: Expected range of hard coral cover.
#'   \item \strong{hcc_growth}: Annual growth rate of hard coral cover.
#'   \item \strong{sc_cover_range}: Expected range of soft coral cover.
#'   \item \strong{sc_growth}: Annual growth rate of soft coral cover.
#' }
#'
#' **Large-scale sampling parameters (`config_lrge`)**
#' \itemize{
#'   \item \strong{n_locs}: Number of reef locations.
#'   \item \strong{n_sites}: Number of sites per reef.
#'   \item \strong{seed}: Seed for sampling reproducibility.
#' }
#'
#' **Fine-scale sampling parameters (`config_fine`)**
#'
#' Defines site-level, transect-level, and residual variance components for
#' each benthic group.
#'
#' \itemize{
#'   \item \strong{years}: Years included in the fine-scale simulation.
#'   \item \strong{Number_of_transects_per_site}: Transects simulated at each site.
#'   \item \strong{Depths}: Number of depth strata.
#'   \item \strong{hcc_site_sigma}, \strong{hcc_transect_sigma}, \strong{hcc_sigma}:  
#'     Random-effect variances for hard coral.
#'   \item \strong{sc_site_sigma}, \strong{sc_transect_sigma}, \strong{sc_sigma}:  
#'     Random-effect variances for soft coral.
#'   \item \strong{ma_site_sigma}, \strong{ma_transect_sigma}, \strong{ma_sigma}:  
#'     Random-effect variances for macroalgae.
#' }
#'
#' **Point-based sampling parameters (`config_pt`)**
#'
#' Defines the structure of quadrats, frames, and annotation points.
#'
#' \itemize{
#'   \item \strong{Depths}: Number of depth strata.
#'   \item \strong{Depth_effect_multiplier}: Strength of depth-related differences.
#'   \item \strong{Number_of_transects_per_site}: Transects per site.
#'   \item \strong{Number_of_frames_per_transect}: Frames per transect.
#'   \item \strong{Number_of_quadrats_per_transect}: Quadrats per transect.
#'   \item \strong{Points_per_frame}: Annotation points per frame.
#'   \item \strong{Quad_sigma}: Quadrats-level random variation.
#' }
#'
#' @return Invisibly returns a list containing the four configuration lists:
#' \code{config_sp}, \code{config_lrge}, \code{config_fine}, and
#' \code{config_pt}. Each is also assigned to the global environment.
#'
#' @author Murray
#'
#' @export
generateSettings <- function(nreefs, nsites, nyears, dhw_eff, cyc_eff, other_eff){

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
  dhw_weight = dhw_eff,
  cyc_weight = cyc_eff,
  other_weight = other_eff,
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
