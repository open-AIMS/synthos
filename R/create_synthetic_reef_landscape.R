#' Create Synthetic Reef Landscape
#'
#' Generates a synthetic reef landscape from a spatial grid and configuration parameters.
#' This includes generating synthetic fields, patches, and reefs, SPDE projections,
#' baseline and synthetic benthos layers (hard coral, soft coral, macroalgae),
#' and optionally reef-level disturbance layers (CYC, DHW, OTHER). Reef coordinates are
#' preserved in the output data frame.
#'
#' @title Create Synthetic Reef Landscape with Benthos and Disturbances
#'
#' @param spatial_grid An sf object representing the spatial grid over which to simulate reefs.
#' @param config A list of configuration parameters controlling the simulation, including:
#'               - `seed`: random seed
#'               - `model`, `psill`, `range`, `nugget`: variogram/SPDE parameters
#'               - `patch_threshold`, `reef_width`: reef patch parameters
#'               - `years`: numeric vector of years to simulate
#'               - `dhw_weight`, `cyc_weight`, `other_weight`: weights for disturbance layers
#'               - `hcc_cover_range`, `hcc_growth`: baseline range and growth rate of hard coral
#'               - `sc_cover_range`, `sc_growth`: baseline range and growth rate of soft coral
#' @param include_disturbances Logical; if TRUE, reef-level disturbances (CYC, DHW, OTHER) are included in the output.
#' @param verbose Logical; if TRUE, prints progress messages during simulation.
#'
#' @return A data.frame containing reef-level benthos data (`HCC`, `SC`, `MA`) and, if `include_disturbances = TRUE`,
#'         disturbance values (`CYC`, `DHW`, `OTHER`), along with reef coordinates.
#'
#' @author Murray
#' @export
create_synthetic_reef_landscape <- function(spatial_grid, config, include_disturbances = FALSE, verbose = FALSE) {
  
  testthat::expect(
    inherits(spatial_grid, c("sfc")),
    "spatial_grid must be an sfc object"
  )
  testthat::expect_in(
    sort(c(
      "seed",
      "model",
      "psill",
      "range",
      "nugget",
      "patch_threshold",
      "reef_width",
      "years",
      "alpha",
      "kappa",
      "variance",
      "dhw_weight",
      "cyc_weight",
      "other_weight",
      "hcc_cover_range",
      "hcc_growth",
      "sc_cover_range",
      "sc_growth"
    )),
    sort(names(config))
  )
  if (verbose) cat("Generating synthetic field\n")
  simulated_field <- synthos::generate_field(spatial_grid, config)
  if (verbose) cat("Generating synthetic patches\n")
  simulated_patches <- synthos::generate_patches(simulated_field, config)
  if (verbose) cat("Generating synthetic reefs\n")
  simulated_reefs <- synthos::generate_reefs(simulated_patches, config)
  if (verbose) cat("Generating SPDE\n")
  matern_projection <- synthos::create_spde(spatial_grid, config)
  if (verbose) cat("Generating DHW layer\n")
  dhw <- synthos::disturbance_dhw(spatial_grid, matern_projection, config)
  if (verbose) cat("Generating cyclone layer\n")
  cyc <- synthos::disturbance_cyc(spatial_grid, matern_projection, config)
  if (verbose) cat("Generating other disturbance layer\n")
  other <- synthos::disturbance_other(spatial_grid, matern_projection, config)
  if (verbose) cat("Combine all effect layers\n")
  all_disturbance_effects <- synthos::disturbance_all(
    spatial_grid,
    dhw_effects = dhw$dhw_effects,
    cyc_effects = cyc$cyc_effects,
    other_effects = other$other_effects,
    matern_projection,
    config) 
  if (verbose) cat("Generate baseline hard coral cover\n")
  baseline_hcc <- synthos::baseline_hard_coral_cover(spatial_grid, matern_projection,
    cover_range = config$hcc_cover_range)
  if (verbose) cat("Generate synthetic hard coral cover\n")
  field_hcc <- synthos::synthetic_field_hcc(spatial_grid, all_disturbance_effects$all_effects_df,
    baseline_hcc$baseline_sample_hcc, matern_projection, config)
  if (verbose) cat("Generate baseline soft coral cover\n")
  baseline_sc <- synthos::baseline_soft_coral_cover(spatial_grid, matern_projection,
    cover_range = config$sc_cover_range, config)
  if (verbose) cat("Generate synthetic soft coral cover\n")
  field_sc <- synthos::synthetic_field_sc(spatial_grid, all_disturbance_effects$all_effects_df,
    baseline_sc$baseline_sample_sc, matern_projection, config)
  if (verbose) cat("Pointify polygons\n")
  reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
  if (verbose) cat("Calculate reef hard coral cover\n")
  reefs_hcc <- synthos::calculate_reef_hcc(
    spatial_grid,
    matern_projection,
    field_hcc$all_effects_hcc,
    reefs$data_reefs_df,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf,
    config
  )
  if (verbose) cat("Calculate reef soft coral cover\n")
  reefs_sc <- synthos::calculate_reef_sc(
    spatial_grid,
    matern_projection,
    field_sc$all_effects_sc,
    reefs$data_reefs_df,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf,
    config
  )
  if (verbose) cat("Calculate reef macroalgae cover\n")
  reefs_ma <- synthos::calculate_reef_ma(
    reefs_hcc$data_reefs_hcc,
    reefs_sc$data_reefs_sc,
    reefs$data_reefs_sf,
    simulated_reefs$simulated_reefs_poly_sf
  )
  if (verbose) cat("Combine reef-level benthos data\n")
  benthos_reefs_pts <- synthos::combine_reef_benthos(
    reefs_hcc$data_reefs_pts_hcc_sf,
    reefs_sc$data_reefs_pts_sc_sf,
    reefs_ma$data_reefs_pts_ma_sf
  )
  if (include_disturbances) {
    if (verbose) cat("Calculate reef disturbances\n")
    reefs_cyc <- synthos::calculate_reef_disturbances(
      spatial_grid,
      matern_projection,
      cyc$cyc_effects,
      reefs$data_reefs_df,
      reefs$data_reefs_sf,
      simulated_reefs$simulated_reefs_poly_sf,
      config
    )
    reefs_dhw <- synthos::calculate_reef_disturbances(
      spatial_grid,
      matern_projection,
      dhw$dhw_effects,
      reefs$data_reefs_df,
      reefs$data_reefs_sf,
      simulated_reefs$simulated_reefs_poly_sf,
      config
    )
    reefs_other <- synthos::calculate_reef_disturbances(
      spatial_grid,
      matern_projection,
      other$other_effects,
      reefs$data_reefs_df,
      reefs$data_reefs_sf,
      simulated_reefs$simulated_reefs_poly_sf,
      config
    )
    benthos_reefs_pts <- synthos::combine_reef_disturbances(
      benthos_reefs_pts,
      reefs_cyc$data_reefs_pts_disturb_sf,
      reefs_dhw$data_reefs_pts_disturb_sf,
      reefs_other$data_reefs_pts_disturb_sf
    )
  }
  return(benthos_reefs_pts)
}


