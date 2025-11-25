#' Generate Simulation Settings for Synthetic Reef Landscapes
#'
#' Creates and saves all configuration lists required to simulate synthetic
#' reef landscapes, including spatial, large-scale, fine-scale, and
#' point-based sampling parameters. Also assigns each configuration to the
#' global environment for downstream modelling functions.
#'
#' @title Generate Settings for Synthetic Reef Simulations
#'
#' @param nreefs Integer. Number of reef locations to simulate.
#' @param nyears Integer. Number of years to simulate.
#'
#' @return Invisibly returns a list containing:
#'   \itemize{
#'     \item `config_sp` – spatio-temporal model configuration
#'     \item `config_lrge` – large-scale sampling configuration
#'     \item `config_fine` – fine-scale sampling configuration
#'     \item `config_pt` – point-based sampling configuration
#'   }
#'   Each configuration list is also assigned to the global environment.
#'
#' @details
#' The function prepares:
#' \itemize{
#'   \item Spatio-temporal parameters for the synthetic field and disturbances
#'   \item Large-scale sampling parameters (reef-level)
#'   \item Fine-scale sampling parameters (site and transect-level)
#'   \item Point-based sampling parameters (photo-quadrats)
#' }
#'  
#' All configurations are saved to disk as `lists_of_parameters.RData`.
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
  dhw_weight = 0.6,
  cyc_weight = 0.39,
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
  Depths = 1,
#  Number_of_frames_per_transect = 100,
#  Points_per_frame = 5,
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
#save(config_list, file = paste0(title_of_run,"/lists_of_parameters.RData"))

## Type of monitoring surveys
surveys <- "fixed"
assign("surveys", surveys, envir = .GlobalEnv)
}
