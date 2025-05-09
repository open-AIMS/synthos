##' Generate simulated patches
##'
##' Generate simulated patches from a simulated field.  In this context
##' a patch is defined as an isolated  area of a field.  This could
##' represent a stand of a specific species of tree, or a coral reef as
##' just two examples.  The patches are generated by filtering the field
##' to only include areas that are above a certain threshold.
##' @title Generate simulated patches 
##' @param simulated_field_sf
##' A simulated field as an sf object as generated by the function
##' \code{generate_field}
##' @param config 
##' A list that contains the parameters for filtering the field to patches.
##' The list should contain the following parameters:
##' - patch_threshold: a numeric representing the threshold below which
##'   the field is masked away to leave only the patches
##' @return an sf object that represents the simulated patches 
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
##'  patch_threshold = 1.75
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
##' simulated_patches |> ggplot() + geom_sf(fill="lightblue")
##' @author Murray
##' @export
generate_patches <- function(simulated_field_sf, config) {
  testthat::expect(
    inherits(simulated_field_sf, c("sf")),
    "simulated_field_sf must be an sf object"
  )
  testthat::expect_match(
    names(simulated_field_sf),
    "sim[0-9]*|geometry"
  )
  testthat::expect_in(
    sort(c("patch_threshold")),
    sort(names(config))
  )
  ## create a thresholded field
  sf_use_s2(FALSE)
  simulated_patches_sf <- simulated_field_sf |>
    dplyr::rename(Y = sim1) |> 
    sf::st_as_sf(coords = c("Longitude", "Latitude")) |>
    dplyr::filter(Y > config$patch_threshold) |>
    sf::st_buffer(0.05, endCapStyle = "SQUARE") |>
    sf::st_cast("POLYGON") |>
    sf::st_union() |>
    ## st_set_crs(4326)
    suppressMessages() |>
    suppressWarnings()
  sf_use_s2(TRUE)
  return(simulated_patches_sf)
}
