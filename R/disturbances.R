##' Degree heating weeks disturbance layer
##'
##' This function generates a disturbance layer for degree heating weeks.
##' This layer is created by combining a temporal trend with a spatially
##' varying random field.
##'
##' The temporal trend is a linear trend with a sinusoidal component.
##' \deqn{
##' cyear = year - 1
##' }
##' \deqn{
##' dhw_i ~ Beta(\alpha_i, 1)
##' }
##' \deqn{
##' log(\frac{\alpha_i}{1 - \alpha_i}) = 0.2 + cYear_i + sin(cYear_i)
##' }
##' \deqn{
##' dhw_i = 5 * \frac{dhw_i - min(dhw_i)}{max(dhw_i) - min(dhw_i)}
##' }
##' The spatially varying random field is generated
##' by a random walk with a beta distributed correlation coefficient.
##' The disturbance layer is then projected onto the spatial grid.
##'
##' @title Degree heating weeks disturbance layer
##' @param spatial_grid
##' A sfc POINT object representing the full spatial grid/
##' @param spde
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @param config
##' A list that should contains the following parameters:
##' - years: A vector of years to generate the disturbance layer for
##' - seed: A seed for the random number generator
##' @return
##' A list containing the following elements:
##' - dhw_temporal: A data frame containing the temporal trend in DHW
##' - dhw_effects: A matrix containing the spatially varying random field
##' - dhw_pts_sample: A matrix containing the disturbance layer projected onto the spatial grid
##' - dhw_pts_effects_df: A data frame containing the disturbance layer projected onto the spatial grid
##' @examples
##' library(sf)
##' library(INLA)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
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
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' dhw$dhw_temporal |>
##'   ggplot(aes(y = Y, x = Year)) +
##'   geom_line() +
##'   theme_bw(base_size = 7)
##'
##' ggplot(dhw$dhw_pts_effects_df, aes(y = Latitude, x = Longitude)) +
##'   geom_tile(aes(fill = Value)) +
##'   facet_wrap(~Year, nrow = 2) +
##'   scale_fill_gradientn(colors = rev(heat.colors(10))) +
##'   coord_sf(crs = 4236) +
##'   theme_bw(base_size = 12) +
##'   theme(axis.title = element_blank())
##' @author Murray
##' @export
disturbance_dhw <- function(spatial_grid, spde, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  ## Overall temporal trend in DHW
  set.seed(config$seed)
  dhw_temporal <- data.frame(Year = config$years) |>
    dplyr::mutate(
      cYear = Year - 1, # as.vector(scale(Year, scale=FALSE)),
      Y = 0.2 * cYear + sin(cYear),
      Y = Y * rbeta(length(config$years), Y, 1),
      Y = scales::rescale(Y - min(Y), to = c(0, 5))
    )
  ## Propagate this temporal trend across a random rield with a time
  ## varying autocorrelation coefficient drawn from a beta distribution
  ## with shape parameters of 0.2 and 1
  set.seed(config$seed)
  dhw_sample <- INLA::inla.qsample(length(config$years),
    spde$Q,
    seed = config$seed,
    constr = spde$spde$f$extraconstr
  ) |>
    suppressMessages() |>
    suppressWarnings()

  rho <- rep(0.7, length(config$years))
  rho <- rbeta(length(config$years), 0.2, 1)
  x <- dhw_sample
  for (j in 2:length(config$years)) {
    x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * dhw_sample[, j]
  }
  x <- sweep(x, 2, dhw_temporal$Y, FUN = "+")
  dhw_effects <- scales::rescale(x, to = c(0, 1))
  dhw_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    dhw_effects
  )

  dhw_pts_effects_df <- dhw_pts_sample %>%
    as.matrix() %>%
    as.data.frame() %>%
    cbind(spatial_grid_pts_df) %>%
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) %>%
    dplyr::mutate(Year = config$years[as.numeric(Year)])
    ## dplyr::mutate(Year = as.numeric(Year))
  ## Value=scales::rescale(Value, to=c(0,1)))
  list(
    dhw_temporal = dhw_temporal,
    dhw_effects = dhw_effects,
    dhw_pts_sample = dhw_pts_sample,
    dhw_pts_effects_df = dhw_pts_effects_df
  )
}

##' Cyclones disturbance layer
##'
##' This function generates a disturbance layer for cyclones.
##' This layer is created by combining a temporal trend with a spatially
##' varying random field.
##'
##' For each year, the probability of a cyclone occurrence is calculated
##' for somewhere within the spatial domain.  The cyclone intensity and
##' a sine wave path for the cyclone to follow is then calculated.
##'
##' @title Cyclones disturbance layer
##' @param spatial_grid
##' A sfc POINT object representing the full spatial grid/
##' @param spde
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @param config
##' A list that should contains the following parameters:
##' - years: A vector of years to generate the disturbance layer for
##' - seed: A seed for the random number generator
##' @return
##' A list containing the following elements:
##' - cyc_temporal: A data frame containing the temporal trend in CYC
##' - cyc_effects: A matrix containing the spatially varying random field
##' - cyc_pts_sample: A matrix containing the disturbance layer projected onto the spatial grid
##' - cyc_pts_effects_df: A data frame containing the disturbance layer projected onto the spatial grid
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
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
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' cyc$cyc_effects_df |>
##'    group_by(Year) |>
##'     summarise(
##'       Mean = mean(Value),
##'       Median = median(Value)
##'     ) |>
##'     ggplot(aes(x = Year)) +
##'     geom_line(aes(y = Mean), color = "blue") +
##'     geom_line(aes(y = Median), color = "red") +
##'     theme_bw(base_size = 7)
##' 
##'  ggplot(cyc$cyc_pts_effects, aes(y = Latitude, x = Longitude)) +
##'    geom_tile(aes(fill = Value)) +
##'    facet_wrap(~Year, nrow = 2) +
##'    scale_fill_gradientn(colors = terrain.colors(10)) +
##'    coord_sf(crs = 4236) +
##'    theme_bw(base_size = 12) +
##'    theme(
##'      axis.title = element_blank(),
##'      legend.position = c(0.95, 0.95),
##'      legend.justification = c(1, 1)
##'    )
##' @author Murray
##' @export
disturbance_cyc <- function(spatial_grid, spde, config) {

  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  set.seed(config$seed)
  cyc <- vector("list", length(config$years))

  yrs <- 1 + (config$years - min(config$years))
  ## for (yr in config$years) {
  for (yr in yrs) {
    ## cat(paste("Year:", yr, "\n"))
    cyc_occur <- rbinom(1, 1, prob = min(0.05 * yr^2, 0.6))
    ## cat(paste("Cyclone Occurance:", cyc_occur, "\n"))
    cyc_intensity <- rbeta(1, 2, 1) |> round(2)
    ## cat(paste("Cyclone intensity:", cyc_intensity, "\n"))
    ## cyc_spatial <- spatial_grid_pts_df  |>
    lat_offset <- runif(1, 0, 5)
    cyc_spatial <- spde$mesh$loc[, 1:2] |>
      as.data.frame() |>
      dplyr::select(Longitude = V1, Latitude = V2) |>
      dplyr::mutate(
        clong = as.vector(scale(Longitude, scale = FALSE)),
        clat = as.vector(scale(Latitude, scale = FALSE)),
        Y = lat_offset + runif(1, -1, 1) * clong + runif(1, -1, 1) *
          clat + sin(clat),
        # Y= Y - runif(1,-10,10),
        Y = abs(Y),
        Y = ifelse(Y > cyc_intensity, cyc_intensity, Y),
        Y = cyc_intensity - Y,
        Value = Y * cyc_occur
      )
    cyc[[yr]] <- cyc_spatial |>
      dplyr::mutate(Year = yr)
  }
  cyc <- do.call("rbind", cyc)
  cyc_effects_df <- cyc |>
    dplyr::mutate(Value = scales::rescale(Value, to = c(0, 1)))

  cyc_effects <- cyc_effects_df |>
    dplyr::select(-clong, -clat, -Y) |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = Value
    ) |>
    dplyr::select(-Longitude, -Latitude) |>
    as.matrix()

  cyc_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    cyc_effects
  )

  cyc_pts_effects <- cyc_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = config$years[as.numeric(Year)])
  list(
    cyc_effects = cyc_effects,
    cyc_effects_df = cyc_effects_df,
    cyc_pts_sample = cyc_pts_sample,
    cyc_pts_effects = cyc_pts_effects
  )
}


##' Other disturbance layer
##'
##' This function generates a disturbance layer for other effects.
##' This layer is created by combining a temporal trend with a spatially
##' varying random field.
##'
##' Similar to degree heating weeks, other disturbances are generated by
##' defining local effects that are autocorrelated. This is a catchall
##' for all other disturbances including crown of thorns, disease etc.
##'
##' @title Others disturbance layer
##' @param spatial_grid
##' A sfc POINT object representing the full spatial grid/
##' @param spde
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @param config
##' A list that should contains the following parameters:
##' - years: A vector of years to generate the disturbance layer for
##' - seed: A seed for the random number generator
##' @return
##' A list containing the following elements:
##' - other_temporal: A data frame containing the temporal trend in OTHER
##' - other_effects: A matrix containing the spatially varying random field
##' - other_pts_sample: A matrix containing the disturbance layer projected onto the spatial grid
##' - other_pts_effects_df: A data frame containing the disturbance layer projected onto the spatial grid
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
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
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' other$other_pts_effects |>
##'  group_by(Year) |>
##'  summarise(
##'    Mean = mean(Value, na.rm = TRUE),
##'    Median = median(Value, na.rm = TRUE)
##'  ) |>
##'  ggplot(aes(x = Year)) +
##'  geom_line(aes(y = Mean))
##'
##' ggplot(other$other_pts_effects, aes(y = Latitude, x = Longitude)) +
##'  geom_tile(aes(fill = Value)) +
##'  facet_wrap(~Year, nrow = 2) +
##'  scale_fill_gradientn(colors = terrain.colors(10)) +
##'  coord_sf(crs = 4236) +
##'  theme_bw(base_size = 12) +
##'  theme(axis.title = element_blank())
##' @author Murray
##' @export
disturbance_other <- function(spatial_grid, spde, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)
  set.seed(config$seed + 1)
  other_sample <- INLA::inla.qsample(length(config$years),
    spde$Q,
    seed = config$seed + 1,
    constr = spde$spde$f$extraconstr
  ) |>
    suppressMessages() |>
    suppressWarnings()

  rho <- rep(0.7, length(config$years))
  rho <- rbeta(length(config$years), 0.2, 1)
  x <- other_sample
  for (j in 2:length(config$years)) {
    x[, j] <- rho[j] * x[, j - 1] + sqrt(1 - rho[j]^2) * other_sample[, j]
  }
  other_effects <- scales::rescale(x, to = c(0, 1))
  other_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    other_effects
  )

  other_pts_effects <- other_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    ## dplyr::mutate(Year = as.numeric(Year)) # ,
    dplyr::mutate(Year = config$years[as.numeric(Year)])
  ## Value=scales::rescale(Value, to=c(0,1)))

  list(
    other_effects = other_effects,
    other_pts_sample = other_pts_sample,
    other_pts_effects = other_pts_effects
  )
}

##' All disturbance summarise_layers()
##'
##' This function combines all the effects together (weighted).
##' The total effects are compiled (disturbances as well as growth) together.
##' All calculations are on the link scale, so as to simply the
##' calculations of a cumulative sum of effects per pixel.
##' The relative influence (annual decline weighting) of each disturbance
##' is provided via the config.
##' Growth (annual increase) of HCC and SC are also provided in the config.
##' Macroalgae respond differently. Rather than respond directly, macroalgae
##' simply takes up the remaining available space \eqn{MA = Total available space - HCC - SC}
##' @title All disturbance layers
##' @param spatial_grid
##' A sfc POINT object representing the full spatial grid/
##' @param dhw_effects
##' A matrix of dhw effect parameters
##' @param cyc_effects
##' A matrix of cyclone effect parameters
##' @param other_effects
##' A matrix of other effect parameters
##' @param spde
##' A list containing the SPDE mesh, SPDE object, Q matrix and A matrix
##' @param config
##' A list that should contains the following parameters:
##' - years: A vector of years to generate the disturbance layer for
##' - seed: A seed for the random number generator
##' - dhw_weight: the relative influence of DHW as a disturbance
##' - cyc_weight: the relative influence of Cyclones as a disturbance
##' - other_weight: the relative influence of other disturbances as a disturbance
##' - hcc_growth: the annual growth rate of hard coral
##' - sc_growth: the annual growth rate of soft coral
##' @return
##' A list containing the following elements:
##' - other_temporal: A data frame containing the temporal trend in OTHER
##' - other_effects: A matrix containing the spatially varying random field
##' - other_pts_sample: A matrix containing the disturbance layer projected onto the spatial grid
##' - other_pts_effects_df: A data frame containing the disturbance layer projected onto the spatial grid
##' @examples
##' library(sf)
##' library(INLA)
##' library(dplyr)
##' library(ggplot2)
##' config <- list(crs=4326, seed = 123)
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
##' config <- list(alpha = 2, kappa = 1, variance = 1)
##' matern_projection <- create_spde(spatial_grid, config)
##' config <- list(years = 1:12, seed = 123)
##' dhw <- disturbance_dhw(spatial_grid, matern_projection, config)
##' cyc <- disturbance_cyc(spatial_grid, matern_projection, config)
##' other <- disturbance_other(spatial_grid, matern_projection, config)
##' config <- list(
##'   years = 1:12, seed = 123,
##'   dhw_weight = 0.5,
##'   cyc_weight = 0.4,
##'   other_weight = 0.1,
##'   hcc_growth = 0.3,
##'   sc_growth =  0.3
##' )
##' all_disturbance_effects <- disturbance_all(
##'   spatial_grid,
##'   dhw_effects = dhw$dhw_effects,
##'   cyc_effects = cyc$cyc_effects,
##'   other_effects = other$other_effects,
##'   matern_projection,
##'   config) 
##' all_disturbance_effects$all_effects_df |>
##'    group_by(Year) |>
##'    summarise(
##'      Mean = mean(Y_HCC, na.rm = TRUE),
##'      Median = median(Y_HCC, na.rm = TRUE)
##'    ) |>
##'    ggplot(aes(x = as.numeric(as.character(Year)))) +
##'    geom_line(aes(y = Mean, color = "Mean")) +
##'    geom_line(aes(y = Median, color = "Median")) +
##'    scale_x_continuous("Year") +
##'    scale_y_continuous("Effect on HCC") +
##'    theme_bw()
##'
##' all_disturbance_effects$all_effects_df |>
##'   ggplot(aes(y = Latitude, x = Longitude)) +
##'   geom_point(aes(color = Y_HCC)) +
##'   ## geom_tile(aes(fill = Y_HCC)) +
##'   facet_wrap(~Year, nrow = 2) +
##'   scale_color_gradient2("HCC", low = "red", high = "green", mid = "white") +
##'   coord_sf(crs = 4236) +
##'   theme_bw(base_size = 12) +
##'   theme(axis.title = element_blank())
##'
##' ggplot(all_disturbance_effects$disturb_pts_effects, aes(y = Latitude, x = Longitude)) +
##' geom_tile(aes(fill = Value)) +
##' facet_wrap(~Year, nrow = 2) +
##' scale_fill_gradientn(colors = terrain.colors(10)) +
##' coord_sf(crs = 4236) +
##' theme_bw(base_size = 12) +
##' theme(
##'   axis.title = element_blank(),
##'   legend.position = c(0.95, 0.95),
##'   legend.justification = c(1, 1)
##' )
##' @author Murray
##' @export
disturbance_all <- function(spatial_grid, dhw_effects, cyc_effects, other_effects, spde, config) {
  testthat::expect(
    inherits(spatial_grid, c("sfc_POINT")),
    "spatial_grid must be an sfc_POINT object"
  )
  testthat::expect_in(
    sort(c("mesh", "spde", "Q", "A")),
    sort(names(spde))
  )
  testthat::expect(
    inherits(spde$spde, c("inla.spde")),
    "spde$spde must be a inla.spde object"
  )
  testthat::expect_in(
    sort(c("years", "seed", "dhw_weight", "cyc_weight", "other_weight")),
    sort(names(config))
  )

  spatial_grid_pts_df <- spatial_grid_sfc_to_df(spatial_grid)

  disturb_effects <-
    (config$dhw_weight * dhw_effects) +
    (config$cyc_weight * cyc_effects) +
    (config$other_weight * other_effects) |>
    as.data.frame() # |>
  all_effects_df <- spde$mesh$loc[, 1:2] |>
    as.data.frame() |>
    dplyr::rename(Longitude = V1, Latitude = V2) |>
    cbind(disturb_effects) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = "Year",
      names_pattern = "sample:(.*)",
      values_to = "Y"
    ) |>
    dplyr::mutate(Year = factor(Year, levels = sort(unique(as.numeric(
      as.character(Year)
    ))))) |>
    dplyr::group_by(Longitude, Latitude) |>
    dplyr::mutate(
      Growth_HCC = 0.3, ## Add growth onto this
      Growth_SC = 0.3,
      Y_HCC = cumsum(-Y + Growth_HCC), ## cumsum on link scale will accumulate effects
      Y_SC = cumsum(-Y + Growth_SC)
    )
  all_effects <- all_effects_df |>
    tidyr::pivot_wider(
      id_cols = c(Longitude, Latitude),
      names_prefix = "sample:",
      names_from = Year,
      values_from = Y_HCC
    )

  ## Project onto the spatial grid
  disturb_pts_sample <- INLA::inla.mesh.project(spde$mesh,
    loc = as.matrix(spatial_grid_pts_df[, 1:2]),
    all_effects |>
      dplyr::ungroup() |> 
      dplyr::select(-Longitude, -Latitude) |>
      as.matrix()
  )
  disturb_pts_effects <- disturb_pts_sample |>
    as.matrix() |>
    as.data.frame() |>
    cbind(spatial_grid_pts_df) |>
    tidyr::pivot_longer(
      cols = c(-Longitude, -Latitude),
      names_to = c("Year"),
      names_pattern = "sample:(.*)",
      values_to = "Value"
    ) |>
    dplyr::mutate(Year = config$years[as.numeric(Year)])
  list(
    disturb_effects = disturb_effects,
    all_effects_df = all_effects_df,
    all_effects =  all_effects,
    disturb_pts_sample = disturb_pts_sample,
    disturb_pts_effects = disturb_pts_effects
  )
}
