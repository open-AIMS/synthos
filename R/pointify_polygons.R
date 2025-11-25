
##' Convert polygon into points
##'
##' Convert a polygon into points.
##' - rasterize the reefs frame
##' - convert to points (centroids of raster cells)
##' - filter to the values of 1
##' - extract coordinates
##' - convert to data frame
##' @title Pointify polygons 
##' @param reefs_sf 
##' A spatial object of class sf with polygons.
##' @return A list with two elements: data_reefs_sf and data_reefs_df.
##' @examples
##' library(sf)
##' library(gstat)
##' library(ggplot2)
##' config <- list(
##'  seed = 1,
##'  crs = 4326,
##'  model = "Exp",
##'  psill = 1,
##'  range = 15,
##'  nugget = 0,
##'  patch_threshold = 1.75,
##'  reef_width = 0.01
##' )
##' spatial_domain <- st_geometry(
##'   st_multipoint(
##'     x = rbind(
##'       c(0, -10),
##'       c(3, -10),
##'       c(10, -20),
##'       c(1, -21),
##'       c(2, -16),
##'       c(0, -10)
##'     )
##'   )
##' ) |>
##'   st_set_crs(config$crs) |>
##'   st_cast("POLYGON")
##' set.seed(config$seed)
##' spatial_grid <- spatial_domain |>
##'   st_set_crs(NA) |>
##'   st_sample(size = 10000, type = "regular") |>
##'   st_set_crs(config$crs)
##' simulated_field <- generate_field(spatial_grid, config)
##' simulated_patches <- generate_patches(simulated_field, config)
##' simulated_reefs <- generate_reefs(simulated_patches, config)
##' reefs <- pointify_polygons(simulated_reefs$simulated_reefs_sf)
##' @author Murray
##' @export
pointify_polygons <- function(reefs_sf) {
  testthat::expect(
    inherits(reefs_sf, c("sfc")),
    "reefs_sf must be a sfc object"
  )
  data_reefs_sf <- reefs_sf |>
    stars::st_as_stars(dx = 0.01) |>  # rasterize
    sf::st_as_sf(as_points = TRUE) |>
    dplyr::filter(values == 1L)

  data_reefs_df <- data_reefs_sf |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(Longitude = X, Latitude = Y)
  list(data_reefs_sf = data_reefs_sf, data_reefs_df = data_reefs_df)
}