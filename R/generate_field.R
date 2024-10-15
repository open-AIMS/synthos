##' Generate a simulated field
##'
##' Generate a simulated field using a variogram model and an sf object.
##' In this context a field represents a spatially continuous surface
##' that represents some sort of landscape.  This could be a forest or
##' marine system or any other.
##' @title Generate a random field 
##' @param spatial_grid
##' An sf object that represents the spatial grid 
##' @param config
##' A list that contains the parameters for the variogram model.
##' The list should contain the following parameters:
##' - seed: an integer that sets the seed for the random number generator
##' - psill: a numeric value that represents the sill of the variogram model
##' - model: a string that represents the variogram model. The following models are supported: "Sph", "Exp", "Gau", "Lin", "Mat", "Ste", "Pen", "Hug", "Hol", "Cor", "Sphlin", "Sphexp", "Sphgaus", "Sphmat", "Sphste", "Sphpen", "Sphhug", "Sphhol", "Sphcor"
##' - range: a numeric value that represents the range of the variogram model
##' - nugget: a numeric value that represents the nugget of the variogram model
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
##'  nugget = 0
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
##' simulated_field |> ggplot() + geom_sf(aes(colour = sim1))
##' @return an sf object that represents the simulated field
##' @author Murray
##' @export
generate_field <- function(spatial_grid, config) {
  testthat::expect(
    inherits(spatial_grid, c("sf", "sfc")),
    "spatial_grid must be an sf object"
  )
  testthat::expect_in(
    sort(c("seed", "psill", "model", "range", "nugget")),
    sort(names(config))
  )

  set.seed(config$seed)
  ## create a variogram model object that can be used to simulate a
  ## random field
  vgm_model <- gstat::vgm(
    psill = config$psill,
    model = config$model,
    range = config$range,
    nugget = config$nugget
  )
  ## Simulate a random field using gstat with the sf object
  sim <- gstat::gstat(formula = z ~ 1, locations = spatial_grid,
    model = vgm_model,
    ## beta: for an intercept only model, it represents the intercept
    beta = 0,
    ## nmax: the number of nearest observations that should be used for a
    ## kriging prediction or simulation, where nearest is defined in
    ## terms of the space of the spatial locations
    nmax = 20,
    ## dummy: consider these data as a dummy variable
    dummy = TRUE
  )
  ## predict the random field
  simulated_field_sf <- predict(sim, newdata = spatial_grid, nsim = 1)
  return(simulated_field_sf)
}
